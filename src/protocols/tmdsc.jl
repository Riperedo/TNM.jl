# Temperature-Modulated DSC (TMDSC) Protocol

"""
    TMDSCProtocol{T<:Real} <: AbstractThermalProtocol

Sinusoidally modulated temperature program used in Temperature-Modulated DSC (TMDSC):
    T(t) = T_start + q0 * t + A_T * sin(ω * t)
    q(t) = dT/dt = q0 + A_T * ω * cos(ω * t)

# Fields
- `T_start::T`: Starting temperature (K).
- `T_end::T`: Final temperature (K).
- `q0::T`: Underlying linear heating/cooling rate (K/s).
- `amplitude::T`: Temperature modulation amplitude A_T (K).
- `period::T`: Modulation period t_p (seconds).
- `omega::T`: Angular frequency ω = 2π / t_p (rad/s).
- `duration::T`: Total scan duration.
"""
struct TMDSCProtocol{T<:Real} <: AbstractThermalProtocol
    T_start::T
    T_end::T
    q0::T
    amplitude::T
    period::T
    omega::T
    duration::T

    function TMDSCProtocol(T_start::Real, T_end::Real, q0::Real, amplitude::Real, period::Real)
        !iszero(q0) || throw(ArgumentError("Underlying rate q0 cannot be zero for ramp TMDSC"))
        Float64(period) > 0.0 || throw(ArgumentError("Period must be positive, got $period"))
        Float64(amplitude) >= 0.0 || throw(ArgumentError("Amplitude must be non-negative, got $amplitude"))
        
        delta_T = Float64(T_end) - Float64(T_start)
        q0_val = Float64(q0)
        if delta_T * q0_val < 0.0
            q0_val = copysign(q0_val, delta_T)
        end
        dur = abs(delta_T / q0_val)
        omega = 2.0 * π / Float64(period)

        p1, p2, p3, p4, p5, p6, p7 = promote(
            Float64(T_start), Float64(T_end), q0_val,
            Float64(amplitude), Float64(period), omega, dur
        )
        new{typeof(p1)}(p1, p2, p3, p4, p5, p6, p7)
    end
end

# Convenience constructor accepting rate in K/min
function tmdsc_protocol_k_min(T_start::Real, T_end::Real, q0_k_min::Real, amplitude::Real, period_s::Real)
    return TMDSCProtocol(T_start, T_end, q0_k_min / 60.0, amplitude, period_s)
end

duration(p::TMDSCProtocol) = p.duration

function temperature_at(p::TMDSCProtocol, t::Real)
    t_clamped = clamp(t, zero(t), p.duration)
    return p.T_start + p.q0 * t_clamped + p.amplitude * sin(p.omega * t_clamped)
end

function heating_rate_at(p::TMDSCProtocol, t::Real)
    if t < zero(t) || t > p.duration
        return zero(p.q0)
    end
    return p.q0 + p.amplitude * p.omega * cos(p.omega * t)
end
