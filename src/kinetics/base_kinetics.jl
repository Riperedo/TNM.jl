# Base interface for relaxation kinetics

"""
    log_relaxation_time(kinetics::AbstractRelaxationKinetics, T::Real, Tf::Real, units::AbstractUnitSystem=PhysicalUnits())

Computes the natural logarithm of the structural relaxation time ln(τ) at actual temperature `T` and fictive temperature `Tf`.
"""
function log_relaxation_time end

"""
    relaxation_time(kinetics::AbstractRelaxationKinetics, T::Real, Tf::Real, units::AbstractUnitSystem=PhysicalUnits())

Computes the structural relaxation time τ at actual temperature `T` and fictive temperature `Tf`.
"""
function relaxation_time(kinetics::AbstractRelaxationKinetics, T::Real, Tf::Real, units::AbstractUnitSystem=PhysicalUnits())
    return exp(log_relaxation_time(kinetics, T, Tf, units))
end
