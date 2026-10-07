# Linear Ramp and Multi-Segment Thermal Protocols

"""
    LinearRampProtocol{T<:Real} <: AbstractThermalProtocol

Single continuous linear heating or cooling ramp from `T_start` to `T_end` with rate `rate` (K/s or K/min).

# Fields
- `T_start::T`: Initial temperature (K or θ).
- `T_end::T`: Final target temperature (K or θ).
- `rate::T`: Heating/cooling rate dT/dt with proper sign (rate > 0 for heating, rate < 0 for cooling).
- `duration::T`: Total duration of the ramp.
"""
struct LinearRampProtocol{T<:Real} <: AbstractThermalProtocol
    T_start::T
    T_end::T
    rate::T
    duration::T

    function LinearRampProtocol(T_start::Real, T_end::Real, rate::Real)
        !iszero(rate) || throw(ArgumentError("Ramp rate cannot be zero"))
        delta_T = T_end - T_start
        actual_rate = (delta_T * rate < 0) ? copysign(Float64(rate), delta_T) : Float64(rate)
        dur = abs(delta_T / actual_rate)
        p1, p2, p3, p4 = promote(Float64(T_start), Float64(T_end), actual_rate, Float64(dur))
        new{typeof(p1)}(p1, p2, p3, p4)
    end
end

LinearRampProtocol(; T_start::Real, T_end::Real, rate::Real) = 
    LinearRampProtocol(T_start, T_end, rate)

# Convenience constructor allowing rate in K/min
function linear_ramp_k_min(T_start::Real, T_end::Real, rate_k_min::Real)
    rate_k_s = rate_k_min / 60.0
    return LinearRampProtocol(T_start, T_end, rate_k_s)
end

duration(p::LinearRampProtocol) = p.duration

function temperature_at(p::LinearRampProtocol, t::Real)
    t_clamped = clamp(t, zero(t), p.duration)
    return p.T_start + p.rate * t_clamped
end

function heating_rate_at(p::LinearRampProtocol, t::Real)
    if t < zero(t) || t > p.duration
        return zero(p.rate)
    end
    return p.rate
end

"""
    RampSegment{T<:Real}

A single segment in a multi-ramp protocol.
"""
struct RampSegment{T<:Real}
    T_target::T
    rate::T # dT/dt with sign
    duration::T
end

"""
    MultiRampProtocol{T<:Real} <: AbstractThermalProtocol

A piecewise linear sequence of temperature ramps.

# Fields
- `T_start::T`: Initial starting temperature.
- `segments::Vector{RampSegment{T}}`: Ordered list of ramp segments.
- `cumulative_times::Vector{T}`: Time boundaries for each segment.
"""
struct MultiRampProtocol{T<:Real} <: AbstractThermalProtocol
    T_start::T
    segments::Vector{RampSegment{T}}
    cumulative_times::Vector{T}

    function MultiRampProtocol(T_start::T, segments::Vector{RampSegment{T}}) where {T<:Real}
        isempty(segments) && throw(ArgumentError("Segments list cannot be empty"))
        cum_times = zeros(T, length(segments) + 1)
        for i in 1:length(segments)
            cum_times[i + 1] = cum_times[i] + segments[i].duration
        end
        new{T}(T_start, segments, cum_times)
    end
end

function MultiRampProtocol(T_start::Real, targets::Vector{<:Real}, rates::Vector{<:Real})
    length(targets) == length(rates) || throw(DimensionMismatch("Targets and rates must have same length"))
    T = promote_type(typeof(float(T_start)), eltype(targets), eltype(rates))
    current_T = T(T_start)
    segs = RampSegment{T}[]
    for i in 1:length(targets)
        target = T(targets[i])
        rate = T(rates[i])
        delta_T = target - current_T
        if iszero(rate)
            # Isothermal hold if target == current_T
            throw(ArgumentError("Rate cannot be zero; use explicit duration for holds"))
        end
        if delta_T * rate < 0
            rate = copysign(rate, delta_T)
        end
        dur = abs(delta_T / rate)
        push!(segs, RampSegment{T}(target, rate, dur))
        current_T = target
    end
    return MultiRampProtocol(T(T_start), segs)
end

duration(p::MultiRampProtocol) = p.cumulative_times[end]

function _find_segment_index(p::MultiRampProtocol, t::Real)
    t_clamped = clamp(t, zero(t), duration(p))
    idx = searchsortedlast(p.cumulative_times, t_clamped)
    return min(idx, length(p.segments))
end

function temperature_at(p::MultiRampProtocol, t::Real)
    if t <= zero(t)
        return p.T_start
    elseif t >= duration(p)
        return p.segments[end].T_target
    end
    idx = _find_segment_index(p, t)
    t_seg_start = p.cumulative_times[idx]
    T_seg_start = (idx == 1) ? p.T_start : p.segments[idx - 1].T_target
    dt = t - t_seg_start
    return T_seg_start + p.segments[idx].rate * dt
end

function heating_rate_at(p::MultiRampProtocol, t::Real)
    if t < zero(t) || t > duration(p)
        return zero(p.T_start)
    end
    idx = _find_segment_index(p, t)
    return p.segments[idx].rate
end
