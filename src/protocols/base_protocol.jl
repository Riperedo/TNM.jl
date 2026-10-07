# Base interface for thermal protocols

"""
    temperature_at(protocol::AbstractThermalProtocol, t::Real)

Returns the sample temperature T(t) at time `t`.
"""
function temperature_at end

"""
    heating_rate_at(protocol::AbstractThermalProtocol, t::Real)

Returns the instantaneous heating/cooling rate q(t) = dT/dt at time `t`.
"""
function heating_rate_at end

"""
    duration(protocol::AbstractThermalProtocol)

Returns the total duration of the thermal protocol in seconds (or dimensionless time in reduced units).
"""
function duration end

"""
    initial_temperature(protocol::AbstractThermalProtocol)

Returns the starting temperature T0 at time t = 0.
"""
initial_temperature(p::AbstractThermalProtocol) = temperature_at(p, zero(Float64))
