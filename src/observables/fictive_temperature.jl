# Observables: Fictive Temperature and Glass Transition Tg Determination

using Statistics

"""
    limiting_fictive_temperature(res::SimulationResult)

Calculates the limiting fictive temperature Tf' (frozen structural state) from a cooling or heating scan
using the Moynihan equal-area construction:
    Tf' = T_high - ∫_{T_low}^{T_high} Cp^N(T) dT
"""
function limiting_fictive_temperature(res::SimulationResult)
    T = res.temperatures
    Tf = res.fictive_temperatures
    N = length(T)
    
    # Find minimum temperature index (turnaround point in cooling-reheating scans)
    min_idx = argmin(T)
    if min_idx == N
        # Pure cooling scan ending at low T: final Tf is the frozen Tf'
        return Tf[end]
    elseif min_idx == 1
        # Pure heating scan: calculate from Moynihan equal-area integral
        cp_n = normalized_heat_capacity(res)
        integral = 0.0
        for i in 2:N
            dT = T[i] - T[i - 1]
            avg_cp = 0.5 * (cp_n[i] + cp_n[i - 1])
            integral += avg_cp * dT
        end
        return T[end] - integral
    else
        # Multi-segment (cooling followed by reheating): Tf' frozen upon cooling is at turnaround
        return Tf[min_idx]
    end
end

"""
    glass_transition_temperature(res::SimulationResult; method::Symbol=:peak)

Determines the glass transition temperature Tg from the normalized heat capacity peak or inflection point.
Supported methods:
- `:peak`: Temperature corresponding to the maximum overshoot in Cp^N.
- `:inflection`: Temperature corresponding to the maximum rate of change (derivative peak).
- `:half_step`: Temperature where Cp^N crosses 0.5.
"""
function glass_transition_temperature(res::SimulationResult; method::Symbol=:peak)
    T_all = res.temperatures
    cp_n_all = normalized_heat_capacity(res)
    
    # If cooling + heating protocol, isolate the heating segment
    min_idx = argmin(T_all)
    idx_range = (min_idx < length(T_all)) ? (min_idx:length(T_all)) : (1:length(T_all))
    
    T = T_all[idx_range]
    cp_n = cp_n_all[idx_range]
    
    if method === :peak
        idx = argmax(cp_n)
        return T[idx]
    elseif method === :half_step
        for i in 2:length(cp_n)
            if (cp_n[i - 1] <= 0.5 && cp_n[i] >= 0.5) || (cp_n[i - 1] >= 0.5 && cp_n[i] <= 0.5)
                # Linear interpolation
                frac = (0.5 - cp_n[i - 1]) / (cp_n[i] - cp_n[i - 1] + 1e-12)
                return T[i - 1] + frac * (T[i] - T[i - 1])
            end
        end
        return T[argmax(cp_n)]
    elseif method === :inflection
        # Derivative of Cp^N
        dCp_dT = zeros(length(T))
        for i in 2:(length(T) - 1)
            dT = T[i + 1] - T[i - 1]
            dCp_dT[i] = (cp_n[i + 1] - cp_n[i - 1]) / (dT + 1e-12)
        end
        return T[argmax(abs.(dCp_dT))]
    else
        throw(ArgumentError("Unknown Tg method: $method. Use :peak, :half_step, or :inflection"))
    end
end
