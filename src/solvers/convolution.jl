# Moynihan-Markovsky Recursive Convolution Solver

"""
    ConvolutionSolver <: AbstractTNMSolver

Numerical solver using the discrete Moynihan-Markovsky recursive history convolution integral:
    T_{f, n} = T_0 + ∑_{j=1}^n ΔT_j * [ 1 - Φ( ∑_{k=j}^n Δt_k / τ_k ) ]

# Fields
- `n_steps::Int`: Default number of discretization intervals across protocol duration (default: 2000).
- `dt::Union{Nothing, Float64}`: Fixed time step if specified.
"""
struct ConvolutionSolver <: AbstractTNMSolver
    n_steps::Int
    dt::Union{Nothing, Float64}

    function ConvolutionSolver(; n_steps::Int=2000, dt=nothing)
        n_steps > 1 || throw(ArgumentError("n_steps must be greater than 1"))
        new(n_steps, dt === nothing ? nothing : Float64(dt))
    end
end

function solve(model::TNMModel, protocol::AbstractThermalProtocol, solver::ConvolutionSolver;
               units::AbstractUnitSystem=PhysicalUnits(), initial_Tf::Union{Nothing, Real}=nothing)
    dur = duration(protocol)
    dur > zero(dur) || throw(ArgumentError("Protocol duration must be positive"))

    # Construct time grid
    times = if solver.dt !== nothing
        collect(range(0.0, dur, step=solver.dt))
    else
        collect(range(0.0, dur, length=solver.n_steps))
    end
    if times[end] < dur
        push!(times, dur)
    end
    N = length(times)

    T_arr = [temperature_at(protocol, t) for t in times]
    q_arr = [heating_rate_at(protocol, t) for t in times]
    Tf_arr = zeros(Float64, N)

    # Initial condition: equilibrium T_f(0) = T(0) unless specified otherwise
    T0 = T_arr[1]
    Tf0 = initial_Tf === nothing ? T0 : Float64(initial_Tf)
    Tf_arr[1] = Tf0

    # Dispatch based on spectrum type for maximum speed
    _solve_convolution_loop!(Tf_arr, T_arr, times, model.kinetics, model.spectrum, units)

    return SimulationResult(times, T_arr, Tf_arr, q_arr, model, protocol, units)
end

function _solve_convolution_loop!(Tf_arr::Vector{Float64}, T_arr::Vector{Float64}, times::Vector{Float64},
                                 kinetics::AbstractRelaxationKinetics, spectrum::KWWSpectrum, units::AbstractUnitSystem)
    N = length(times)
    beta = spectrum.beta
    T0 = T_arr[1]
    
    # z_accum[j] stores the accumulated material time from step j to current step n
    z_accum = zeros(Float64, N)
    
    @inbounds for n in 2:N
        dt_n = times[n] - times[n - 1]
        
        # Calculate relaxation time at Tk, Tf_{k-1}
        tau_n = relaxation_time(kinetics, T_arr[n], Tf_arr[n - 1], units)
        dz_n = dt_n / tau_n
        
        # Update accumulated material times for all previous temperature steps
        for j in 1:(n - 1)
            z_accum[j] += dz_n
        end
        # The step from n-1 to n has accumulated time dz_n
        z_accum[n] = dz_n

        # Compute fictive temperature via Boltzmann superposition sum
        sum_val = 0.0
        for j in 2:n
            delta_T_j = T_arr[j] - T_arr[j - 1]
            # Φ(z) = exp( - z^β )
            phi_val = exp(-(z_accum[j]^beta))
            sum_val += delta_T_j * (1.0 - phi_val)
        end
        Tf_arr[n] = T0 + sum_val
    end
    return nothing
end

function _solve_convolution_loop!(Tf_arr::Vector{Float64}, T_arr::Vector{Float64}, times::Vector{Float64},
                                 kinetics::AbstractRelaxationKinetics, spectrum::PronySpectrum{M, T_spec}, units::AbstractUnitSystem) where {M, T_spec}
    N = length(times)
    g = spectrum.g
    lambdas = spectrum.lambda
    
    # Internal state variables for each mode Tf_i
    Tf_modes = zeros(Float64, M)
    Tf_modes .= Tf_arr[1]
    
    @inbounds for n in 2:N
        dt_n = times[n] - times[n - 1]
        T_curr = T_arr[n]
        Tf_prev = Tf_arr[n - 1]
        
        tau_curr = relaxation_time(kinetics, T_curr, Tf_prev, units)
        dz = dt_n / tau_curr
        
        # Exact linear relaxation step for each Prony mode:
        # Tf_i(n) = Tf_i(n-1) * exp(-dz / λ_i) + T_curr * (1 - exp(-dz / λ_i))
        current_Tf = 0.0
        for i in 1:M
            decay = exp(-dz / lambdas[i])
            Tf_modes[i] = Tf_modes[i] * decay + T_curr * (1.0 - decay)
            current_Tf += g[i] * Tf_modes[i]
        end
        Tf_arr[n] = current_Tf
    end
    return nothing
end
