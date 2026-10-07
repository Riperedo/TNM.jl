module TNM

using LinearAlgebra
using Statistics
using StaticArrays
using DelimitedFiles
using Printf
using Interpolations
using OrdinaryDiffEq
using Random

# 1. Type hierarchy and unit systems
include("types.jl")

# 2. Kinetics representations
include("kinetics/base_kinetics.jl")
include("kinetics/arrhenius.jl")
include("kinetics/vfth.jl")
include("kinetics/adams_gibbs.jl")

# 3. Relaxation spectrum representations
include("spectrum/kww.jl")
include("spectrum/prony.jl")
include("spectrum/prony_fitting.jl")

# 4. Thermal protocols
include("protocols/base_protocol.jl")
include("protocols/ramps.jl")
include("protocols/isothermal.jl")
include("protocols/tmdsc.jl")
include("protocols/custom.jl")

# 5. Numerical solvers
include("solvers/base_solver.jl")
include("solvers/convolution.jl")
include("solvers/ode_solver.jl")

# 6. Physical observables and analysis
include("observables/heat_capacity.jl")
include("observables/fictive_temperature.jl")
include("observables/tmdsc_analysis.jl")

# 7. Data export and Gnuplot scripting
include("export/dat_exporter.jl")

# 8. Loss functions and metaheuristic parameter estimation
include("fitting/objectives.jl")
include("fitting/differential_evolution.jl")
include("fitting/particle_swarm.jl")

# 9. Experimental data ingestion
include("data/experimental_loader.jl")

# Exports
export
    # Unit Systems & Base Types
    AbstractUnitSystem, PhysicalUnits, ReducedUnits,
    AbstractStructuralRelaxationModel, AbstractRelaxationKinetics,
    AbstractRelaxationSpectrum, AbstractThermalProtocol,
    AbstractTNMSolver, AbstractTNMResult, AbstractOptimizationEngine,
    TNMModel, SimulationResult,

    # Kinetics
    ArrheniusKinetics, VFTHKinetics, AdamsGibbsKinetics,
    relaxation_time, log_relaxation_time,

    # Spectrum
    KWWSpectrum, PronySpectrum, relaxation_function, kww_to_prony,

    # Protocols
    LinearRampProtocol, MultiRampProtocol, RampSegment, linear_ramp_k_min,
    IsothermalAgingProtocol, KovacsCrossoverProtocol,
    TMDSCProtocol, tmdsc_protocol_k_min, TabularProtocol,
    temperature_at, heating_rate_at, duration, initial_temperature,

    # Solvers
    ConvolutionSolver, PronyODESolver, solve,

    # Observables
    normalized_heat_capacity, heat_capacity, enthalpy, heat_flow,
    limiting_fictive_temperature, glass_transition_temperature,
    TMDSCResult, analyze_tmdsc,

    # Export & Visualization
    export_dat, generate_gnuplot_script,

    # Parameter Estimation & Fitting
    ExperimentalDSCData, create_tnm_cost_function,
    TopKTracker, DEResult, optimize_differential_evolution,
    PSOResult, optimize_particle_swarm,

    # Experimental Ingestion
    load_experimental_file, load_experimental_directory

end # module TNM
