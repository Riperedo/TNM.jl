# Base solver definitions, model containers, and simulation results

"""
    TNMModel{K<:AbstractRelaxationKinetics, S<:AbstractRelaxationSpectrum} <: AbstractStructuralRelaxationModel

Complete Tool-Narayanaswamy-Moynihan structural relaxation model combining temperature/structure kinetics and relaxation spectrum.

# Fields
- `kinetics::K`: Relaxation kinetics instance (Arrhenius, VFTH, Adams-Gibbs).
- `spectrum::S`: Relaxation spectrum instance (KWW, Prony).
"""
struct TNMModel{K<:AbstractRelaxationKinetics, S<:AbstractRelaxationSpectrum} <: AbstractStructuralRelaxationModel
    kinetics::K
    spectrum::S
end

TNMModel(kinetics::AbstractRelaxationKinetics) = TNMModel(kinetics, KWWSpectrum(1.0))

"""
    SimulationResult{T<:Real, M<:AbstractStructuralRelaxationModel, P<:AbstractThermalProtocol, U<:AbstractUnitSystem} <: AbstractTNMResult

Result of a structural relaxation simulation.

# Fields
- `times::Vector{T}`: Time points (seconds or dimensionless time).
- `temperatures::Vector{T}`: Sample temperatures T(t) (Kelvin or θ).
- `fictive_temperatures::Vector{T}`: Calculated fictive temperatures Tf(t) (Kelvin or θ_f).
- `heating_rates::Vector{T}`: Instantaneous heating/cooling rates q(t).
- `model::M`: Model used for the simulation.
- `protocol::P`: Thermal history protocol applied.
- `units::U`: Unit system used in calculation.
"""
struct SimulationResult{T<:Real, M<:AbstractStructuralRelaxationModel, P<:AbstractThermalProtocol, U<:AbstractUnitSystem} <: AbstractTNMResult
    times::Vector{T}
    temperatures::Vector{T}
    fictive_temperatures::Vector{T}
    heating_rates::Vector{T}
    model::M
    protocol::P
    units::U
end

"""
    solve(model::AbstractStructuralRelaxationModel, protocol::AbstractThermalProtocol, solver::AbstractTNMSolver; kwargs...)

Executes structural relaxation simulation of `model` under `protocol` using numerical `solver`.
"""
function solve end
