# Tabular and Custom Discrete Thermal Protocols

using Interpolations

"""
    TabularProtocol{T<:Real} <: AbstractThermalProtocol

Thermal protocol defined by discrete experimental or user-supplied time-temperature data series.

# Fields
- `times::Vector{T}`: Time points in ascending order.
- `temperatures::Vector{T}`: Corresponding temperatures.
- `itp_T`: Continuous interpolation object for T(t).
"""
struct TabularProtocol{T<:Real} <: AbstractThermalProtocol
    times::Vector{T}
    temperatures::Vector{T}
    itp_T

    function TabularProtocol(times::Vector{T}, temperatures::Vector{T}) where {T<:Real}
        length(times) == length(temperatures) || throw(DimensionMismatch("times and temperatures must have same length"))
        issorted(times) || throw(ArgumentError("times must be strictly sorted in ascending order"))
        length(times) >= 2 || throw(ArgumentError("Must have at least 2 points"))
        
        itp = linear_interpolation(times, temperatures, extrapolation_bc=Flat())
        new{T}(times, temperatures, itp)
    end
end

TabularProtocol(times::AbstractVector{<:Real}, temperatures::AbstractVector{<:Real}) = 
    TabularProtocol(Vector{Float64}(times), Vector{Float64}(temperatures))

duration(p::TabularProtocol) = p.times[end] - p.times[1]

function temperature_at(p::TabularProtocol, t::Real)
    t_shifted = p.times[1] + clamp(t, zero(t), duration(p))
    return p.itp_T(t_shifted)
end

function heating_rate_at(p::TabularProtocol, t::Real)
    t_shifted = p.times[1] + clamp(t, zero(t), duration(p))
    idx = searchsortedlast(p.times, t_shifted)
    idx = clamp(idx, 1, length(p.times) - 1)
    dt = p.times[idx + 1] - p.times[idx]
    dT = p.temperatures[idx + 1] - p.temperatures[idx]
    return dT / dt
end
