# State-Space ODE Solver for Multi-Mode Structural Relaxation

using OrdinaryDiffEq
using StaticArrays

"""
    PronyODESolver <: AbstractTNMSolver

High-performance stiff ODE solver integrating the coupled non-linear state-space equations:
    d(T_{f,i})/dt = - (T_{f,i}(t) - T(t)) / (λ_i * τ(T(t), Tf(t)))
    Tf(t) = ∑_{i=1}^N g_i * T_{f,i}(t)

# Fields
- `modes::Int`: Number of Prony modes used if KWW spectrum is converted (default: 10).
- `alg`: OrdinaryDiffEq algorithm (default: `Rosenbrock23()`).
- `reltol::Float64`: Relative tolerance (default: 1e-6).
- `abstol::Float64`: Absolute tolerance (default: 1e-6).
- `saveat::Union{Nothing, Vector{Float64}, Float64}`: Time points to save solution at.
"""
struct PronyODESolver <: AbstractTNMSolver
    modes::Int
    alg
    reltol::Float64
    abstol::Float64
    saveat

    function PronyODESolver(; modes::Int=10, alg=Rosenbrock23(), reltol::Float64=1e-6, abstol::Float64=1e-6, saveat=nothing)
        modes > 0 || throw(ArgumentError("modes must be positive"))
        new(modes, alg, reltol, abstol, saveat)
    end
end

function solve(model::TNMModel, protocol::AbstractThermalProtocol, solver::PronyODESolver;
               units::AbstractUnitSystem=PhysicalUnits(), initial_Tf::Union{Nothing, Real}=nothing)
    dur = duration(protocol)
    dur > zero(dur) || throw(ArgumentError("Protocol duration must be positive"))

    # If spectrum is KWW, convert to Prony
    prony_spec = if model.spectrum isa KWWSpectrum
        kww_to_prony(model.spectrum; N=solver.modes)
    elseif model.spectrum isa PronySpectrum
        model.spectrum
    else
        throw(ArgumentError("Unsupported spectrum type $(typeof(model.spectrum)) for ODE solver"))
    end

    M = length(prony_spec.g)
    g_weights = Vector{Float64}(prony_spec.g)
    lambda_shifts = Vector{Float64}(prony_spec.lambda)

    T0 = initial_temperature(protocol)
    Tf0 = initial_Tf === nothing ? T0 : Float64(initial_Tf)
    u0 = fill(Tf0, M)

    # ODE vector field
    function tnm_ode!(du, u, p, t)
        T_curr = temperature_at(protocol, t)
        # Compute macroscopic fictive temperature
        Tf_curr = dot(g_weights, u)
        
        # Calculate relaxation time
        tau_curr = relaxation_time(model.kinetics, T_curr, Tf_curr, units)
        
        @inbounds for i in 1:M
            du[i] = -(u[i] - T_curr) / (lambda_shifts[i] * tau_curr)
        end
        return nothing
    end

    tspan = (0.0, dur)
    save_grid = if solver.saveat isa Vector{Float64}
        solver.saveat
    elseif solver.saveat isa Real
        collect(range(0.0, dur, step=Float64(solver.saveat)))
    else
        collect(range(0.0, dur, length=2000))
    end
    if save_grid[end] < dur
        push!(save_grid, dur)
    end

    prob = ODEProblem(tnm_ode!, u0, tspan)
    sol = OrdinaryDiffEq.solve(prob, solver.alg; reltol=solver.reltol, abstol=solver.abstol, saveat=save_grid)

    times = sol.t
    N_out = length(times)
    T_arr = [temperature_at(protocol, t) for t in times]
    q_arr = [heating_rate_at(protocol, t) for t in times]
    Tf_arr = [dot(g_weights, sol.u[i]) for i in 1:N_out]

    return SimulationResult(times, T_arr, Tf_arr, q_arr, model, protocol, units)
end
