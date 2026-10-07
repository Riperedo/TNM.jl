# Observables: Normalized Heat Capacity, Apparent Heat Capacity, and Enthalpy

"""
    normalized_heat_capacity(res::SimulationResult)

Calculates the apparent normalized heat capacity:
    Cp^N(T) = dTf / dT = (dTf/dt) / (dT/dt)

Returns a vector of normalized heat capacity values corresponding to `res.times`.
"""
function normalized_heat_capacity(res::SimulationResult)
    t = res.times
    T = res.temperatures
    Tf = res.fictive_temperatures
    q = res.heating_rates
    N = length(t)

    dTf_dt = zeros(Float64, N)
    
    # 2nd-order central differences with 1st-order endpoints
    if N >= 3
        dTf_dt[1] = (Tf[2] - Tf[1]) / (t[2] - t[1])
        for i in 2:(N - 1)
            dt_prev = t[i] - t[i - 1]
            dt_next = t[i + 1] - t[i]
            # General non-uniform grid central difference
            dTf_dt[i] = ((Tf[i + 1] - Tf[i]) / dt_next * dt_prev + (Tf[i] - Tf[i - 1]) / dt_prev * dt_next) / (dt_prev + dt_next)
        end
        dTf_dt[N] = (Tf[N] - Tf[N - 1]) / (t[N] - t[N - 1])
    elseif N == 2
        rate = (Tf[2] - Tf[1]) / (t[2] - t[1])
        dTf_dt .= rate
    else
        dTf_dt .= 0.0
    end

    cp_n = zeros(Float64, N)
    for i in 1:N
        if abs(q[i]) > 1e-12
            cp_n[i] = dTf_dt[i] / q[i]
        else
            cp_n[i] = 0.0
        end
    end
    return cp_n
end

"""
    heat_capacity(res::SimulationResult; Cpg::Real=1.0, Cpl::Real=2.0)

Calculates the apparent total heat capacity:
    Cp(T) = Cpg(T) + (Cpl(T) - Cpg(T)) * Cp^N(T)
"""
function heat_capacity(res::SimulationResult; Cpg::Real=1.0, Cpl::Real=2.0)
    cp_n = normalized_heat_capacity(res)
    delta_Cp = Cpl - Cpg
    return Cpg .+ delta_Cp .* cp_n
end

"""
    enthalpy(res::SimulationResult; Cpl::Real=2.0, delta_Cp::Real=1.0)

Calculates the specific enthalpy evolution:
    H(T, Tf) = Cpl * T + delta_Cp * (Tf - T)
"""
function enthalpy(res::SimulationResult; Cpl::Real=2.0, delta_Cp::Real=1.0)
    T = res.temperatures
    Tf = res.fictive_temperatures
    return Cpl .* T .+ delta_Cp .* (Tf .- T)
end

"""
    heat_flow(res::SimulationResult; Cpg::Real=1.0, delta_Cp::Real=1.0)

Calculates the instantaneous specific heat flow rate (dH/dt):
    HF(t) = Cpg * q(t) + delta_Cp * (dTf/dt)
"""
function heat_flow(res::SimulationResult; Cpg::Real=1.0, delta_Cp::Real=1.0)
    t = res.times
    Tf = res.fictive_temperatures
    q = res.heating_rates
    N = length(t)

    dTf_dt = zeros(Float64, N)
    if N >= 2
        dTf_dt[1] = (Tf[2] - Tf[1]) / (t[2] - t[1])
        for i in 2:(N - 1)
            dTf_dt[i] = (Tf[i + 1] - Tf[i - 1]) / (t[i + 1] - t[i - 1])
        end
        dTf_dt[N] = (Tf[N] - Tf[N - 1]) / (t[N] - t[N - 1])
    end

    return Cpg .* q .+ delta_Cp .* dTf_dt
end
