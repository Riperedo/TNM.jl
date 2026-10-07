using Test
using TNM

@testset "TNM.jl Comprehensive Test Suite" begin
    include("test_units_and_kinetics.jl")
    include("test_spectrum.jl")
    include("test_protocols.jl")
    include("test_solvers_and_observables.jl")
    include("test_export.jl")
    include("test_fitting.jl")
    include("test_data_loader.jl")
end
