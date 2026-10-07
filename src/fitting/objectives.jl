# Objective Functions and Loss Formulations for TNM Parameter Estimation

"""
    ExperimentalDSCData

Container for experimental DSC curve data.

# Fields
- `temperatures::Vector{Float64}`: Experimental temperature points.
- `heat_capacities::Vector{Float64}`: Measured normalized heat capacity Cp^N (or apparent Cp).
- `protocol::AbstractThermalProtocol`: Thermal protocol applied during measurement.
- `is_normalized::Bool`: Whether the data is normalized Cp^N or total Cp.
"""
struct ExperimentalDSCData
    temperatures::Vector{Float64}
    heat_capacities::Vector{Float64}
    protocol::AbstractThermalProtocol
    is_normalized::Bool

    function ExperimentalDSCData(temperatures::Vector{Float64}, heat_capacities::Vector{Float64},
                                 protocol::AbstractThermalProtocol; is_normalized::Bool=true)
        length(temperatures) == length(heat_capacities) || throw(DimensionMismatch("Temperatures and heat capacities must have same length"))
        new(temperatures, heat_capacities, protocol, is_normalized)
    end
end

"""
    create_tnm_cost_function(datasets::Vector{ExperimentalDSCData}; 
                             model_type::Symbol=:arrhenius,
                             solver::AbstractTNMSolver=ConvolutionSolver(n_steps=1000),
                             penalty_weight::Float64=1e6)

Creates a loss function `loss(params::AbstractVector{<:Real})` evaluating the Mean Squared Error (MSE)
between modeled heat capacity curves and experimental datasets.

Param vector interpretation:
- `:arrhenius`: `[ln_A, dh_R, x, beta]`
- `:vfth`: `[ln_A, B, T_inf, x, beta]`
"""
function create_tnm_cost_function(datasets::Vector{ExperimentalDSCData}; 
                                 model_type::Symbol=:arrhenius,
                                 solver::AbstractTNMSolver=ConvolutionSolver(n_steps=1000),
                                 penalty_weight::Float64=1e6)
    return function(params::AbstractVector{<:Real})
        # Check boundary bounds and penalties
        if model_type === :arrhenius
            length(params) == 4 || throw(DimensionMismatch("Arrhenius model expects 4 parameters [ln_A, dh_R, x, beta]"))
            ln_A, dh_R, x, beta = params[1], params[2], params[3], params[4]
            if dh_R <= 0 || x <= 0 || x > 1 || beta <= 0 || beta > 1
                return penalty_weight * (1.0 + abs(min(0.0, dh_R)) + abs(min(0.0, x)) + max(0.0, x - 1.0) + abs(min(0.0, beta)) + max(0.0, beta - 1.0))
            end
            kinetics = ArrheniusKinetics(ln_A, dh_R, x)
            spectrum = KWWSpectrum(beta)
            model = TNMModel(kinetics, spectrum)
        elseif model_type === :vfth
            length(params) == 5 || throw(DimensionMismatch("VFTH model expects 5 parameters [ln_A, B, T_inf, x, beta]"))
            ln_A, B, T_inf, x, beta = params[1], params[2], params[3], params[4], params[5]
            if B <= 0 || x <= 0 || x > 1 || beta <= 0 || beta > 1
                return penalty_weight
            end
            kinetics = VFTHKinetics(ln_A, B, T_inf, x)
            spectrum = KWWSpectrum(beta)
            model = TNMModel(kinetics, spectrum)
        else
            throw(ArgumentError("Unknown model_type: $model_type. Use :arrhenius or :vfth"))
        end

        total_loss = 0.0
        total_pts = 0

        for data in datasets
            try
                res = solve(model, data.protocol, solver)
                sim_cp_n = normalized_heat_capacity(res)
                sim_T = res.temperatures

                # Interpolate simulation results onto experimental temperature grid
                # Note: Find nearest or linear interpolation
                for i in 1:length(data.temperatures)
                    T_exp = data.temperatures[i]
                    cp_exp = data.heat_capacities[i]

                    # Find closest simulation temperature index
                    idx = searchsortedlast(sim_T, T_exp)
                    idx = clamp(idx, 1, length(sim_T) - 1)
                    
                    dT = sim_T[idx + 1] - sim_T[idx]
                    frac = abs(dT) > 1e-9 ? (T_exp - sim_T[idx]) / dT : 0.0
                    cp_sim = sim_cp_n[idx] + frac * (sim_cp_n[idx + 1] - sim_cp_n[idx])

                    diff = cp_sim - cp_exp
                    total_loss += diff * diff
                    total_pts += 1
                end
            catch e
                return penalty_weight
            end
        end

        return total_pts > 0 ? (total_loss / total_pts) : penalty_weight
    end
end
