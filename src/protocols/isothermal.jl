# Isothermal Aging and Recovery Protocols

"""
    IsothermalAgingProtocol{T<:Real} <: AbstractThermalProtocol

Full DSC thermal history for aging / structural recovery experiments:
1. Cooling from `T_init` to `T_age` at cooling rate `q_cool` (with q_cool < 0).
2. Isothermal aging hold at `T_age` for time `t_age` (dT/dt = 0).
3. Subsequent reheating from `T_age` to `T_final` at heating rate `q_heat` (with q_heat > 0).

# Fields
- `T_init::T`: Equilibrium starting temperature.
- `T_age::T`: Annealing / aging temperature.
- `T_final::T`: Final temperature after reheating scan.
- `q_cool::T`: Cooling rate (K/s).
- `t_age::T`: Duration of isothermal hold (seconds).
- `q_heat::T`: Reheating rate (K/s).
- `t1::T`: End time of cooling step.
- `t2::T`: End time of isothermal hold step.
- `t3::T`: Total duration (end time of heating step).
"""
struct IsothermalAgingProtocol{T<:Real} <: AbstractThermalProtocol
    T_init::T
    T_age::T
    T_final::T
    q_cool::T
    t_age::T
    q_heat::T
    t1::T
    t2::T
    t3::T

    function IsothermalAgingProtocol(T_init::Real, T_age::Real, T_final::Real,
                                     q_cool::Real, t_age::Real, q_heat::Real)
        Float64(T_init) > Float64(T_age) || throw(ArgumentError("T_init must be greater than T_age for cooling, got T_init=$T_init, T_age=$T_age"))
        Float64(T_final) > Float64(T_age) || throw(ArgumentError("T_final must be greater than T_age for heating, got T_final=$T_final, T_age=$T_age"))
        Float64(t_age) >= 0.0 || throw(ArgumentError("t_age must be non-negative, got $t_age"))
        
        # Ensure correct signs
        q_c = -abs(Float64(q_cool))
        q_h = abs(Float64(q_heat))

        t_cool = abs((Float64(T_age) - Float64(T_init)) / q_c)
        t_hold = Float64(t_age)
        t_heat = abs((Float64(T_final) - Float64(T_age)) / q_h)

        t1 = t_cool
        t2 = t1 + t_hold
        t3 = t2 + t_heat

        p1, p2, p3, p4, p5, p6, p7, p8, p9 = promote(
            Float64(T_init), Float64(T_age), Float64(T_final),
            q_c, t_hold, q_h, t1, t2, t3
        )
        new{typeof(p1)}(p1, p2, p3, p4, p5, p6, p7, p8, p9)
    end
end

duration(p::IsothermalAgingProtocol) = p.t3

function temperature_at(p::IsothermalAgingProtocol, t::Real)
    if t <= zero(t)
        return p.T_init
    elseif t <= p.t1
        # Cooling segment
        return p.T_init + p.q_cool * t
    elseif t <= p.t2
        # Isothermal hold segment
        return p.T_age
    elseif t <= p.t3
        # Heating segment
        dt = t - p.t2
        return p.T_age + p.q_heat * dt
    else
        return p.T_final
    end
end

function heating_rate_at(p::IsothermalAgingProtocol, t::Real)
    if t < zero(t) || t > p.t3
        return zero(p.T_init)
    elseif t <= p.t1
        return p.q_cool
    elseif t <= p.t2
        return zero(p.T_init)
    else
        return p.q_heat
    end
end

"""
    KovacsCrossoverProtocol{T<:Real} <: AbstractThermalProtocol

Protocol for investigating the Kovacs memory/crossover effect:
1. Fast quench from T0 to T_low (annealing temperature).
2. Anneal at T_low for time t_anneal.
3. Rapid jump to intermediate temperature T_target, then observe isothermal relaxation.
"""
struct KovacsCrossoverProtocol{T<:Real} <: AbstractThermalProtocol
    T0::T
    T_low::T
    T_target::T
    t_anneal::T
    t_observe::T
    q_quench::T
    t1::T
    t2::T
    t3::T
    t4::T

    function KovacsCrossoverProtocol(T0::Real, T_low::Real, T_target::Real,
                                     t_anneal::Real, t_observe::Real;
                                     q_quench::Real=100.0)
        t_quench1 = abs(Float64(T0) - Float64(T_low)) / Float64(q_quench)
        t1 = t_quench1
        t2 = t1 + Float64(t_anneal)
        t_jump = abs(Float64(T_target) - Float64(T_low)) / Float64(q_quench)
        t3 = t2 + t_jump
        t4 = t3 + Float64(t_observe)
        
        p1, p2, p3, p4, p5, p6, p7, p8, p9, p10 = promote(
            Float64(T0), Float64(T_low), Float64(T_target),
            Float64(t_anneal), Float64(t_observe), Float64(q_quench),
            t1, t2, t3, t4
        )
        new{typeof(p1)}(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10)
    end
end

duration(p::KovacsCrossoverProtocol) = p.t4

function temperature_at(p::KovacsCrossoverProtocol, t::Real)
    if t <= zero(t)
        return p.T0
    elseif t <= p.t1
        dt = t
        return p.T0 - p.q_quench * dt
    elseif t <= p.t2
        return p.T_low
    elseif t <= p.t3
        dt = t - p.t2
        return p.T_low + p.q_quench * dt
    else
        return p.T_target
    end
end

function heating_rate_at(p::KovacsCrossoverProtocol, t::Real)
    if t <= zero(t) || t > p.t4
        return zero(p.T0)
    elseif t <= p.t1
        return -p.q_quench
    elseif t <= p.t2
        return zero(p.T0)
    elseif t <= p.t3
        return p.q_quench
    else
        return zero(p.T0)
    end
end
