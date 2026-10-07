# Abstract type definitions and unit system representations for TNM.jl

"""
    AbstractUnitSystem

Abstract supertype for unit systems used in structural relaxation calculations.
"""
abstract type AbstractUnitSystem end

"""
    PhysicalUnits <: AbstractUnitSystem

Represents calculations performed using physical SI-compatible units (Temperature in Kelvin, Time in seconds, Rates in K/s or K/min).
"""
struct PhysicalUnits <: AbstractUnitSystem end

"""
    ReducedUnits{T<:Real} <: AbstractUnitSystem

Represents calculations performed using dimensionless reduced variables:
- θ = T / T_ref (dimensionless temperature)
- θ_f = T_f / T_ref (dimensionless fictive temperature)
- t̃ = t / tau_ref (dimensionless time)
- ε = Δh* / (R * T_ref) (dimensionless activation energy)
- q̃ = q * (tau_ref / T_ref) (dimensionless cooling/heating rate)
"""
struct ReducedUnits{T<:Real} <: AbstractUnitSystem
    T_ref::T
    tau_ref::T

    function ReducedUnits(T_ref::Real, tau_ref::Real=1.0)
        T_ref > 0 || throw(ArgumentError("T_ref must be strictly positive, got $T_ref"))
        tau_ref > 0 || throw(ArgumentError("tau_ref must be strictly positive, got $tau_ref"))
        T1, T2 = promote(Float64(T_ref), Float64(tau_ref))
        new{typeof(T1)}(T1, T2)
    end
end

"""
    AbstractStructuralRelaxationModel

Abstract supertype for all glass transition and structural relaxation models.
"""
abstract type AbstractStructuralRelaxationModel end

"""
    AbstractRelaxationKinetics

Abstract supertype for relaxation time temperature/structure dependencies τ(T, Tf).
"""
abstract type AbstractRelaxationKinetics end

"""
    AbstractRelaxationSpectrum

Abstract supertype for relaxation response spectra (e.g. KWW, Prony distribution).
"""
abstract type AbstractRelaxationSpectrum end

"""
    AbstractThermalProtocol

Abstract supertype for temperature-time histories T(t).
"""
abstract type AbstractThermalProtocol end

"""
    AbstractTNMSolver

Abstract supertype for numerical integration solvers for the TNM equations.
"""
abstract type AbstractTNMSolver end

"""
    AbstractTNMResult

Abstract supertype for simulation output results.
"""
abstract type AbstractTNMResult end

"""
    AbstractOptimizationEngine

Abstract supertype for parameter estimation algorithms (e.g., DE, PSO).
"""
abstract type AbstractOptimizationEngine end
