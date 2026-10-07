using Test
using TNM
using StaticArrays
using LinearAlgebra

@testset "Relaxation Spectrum" begin
    @testset "KWW Spectrum" begin
        kww = KWWSpectrum(0.5)
        @test kww.beta == 0.5

        @test relaxation_function(kww, 0.0) == 1.0
        @test isapprox(relaxation_function(kww, 1.0), exp(-1.0^0.5), atol=1e-12)
        @test isapprox(relaxation_function(kww, 4.0), exp(-2.0), atol=1e-12)

        # Debye special case beta = 1.0
        debye = KWWSpectrum(1.0)
        @test isapprox(relaxation_function(debye, 2.5), exp(-2.5), atol=1e-12)

        # Invalid bounds
        @test_throws ArgumentError KWWSpectrum(0.0)
        @test_throws ArgumentError KWWSpectrum(1.5)
        @test_throws ArgumentError relaxation_function(kww, -1.0)
    end

    @testset "Prony Spectrum" begin
        g = SVector{3, Float64}(0.2, 0.5, 0.3)
        lambda = SVector{3, Float64}(0.1, 1.0, 10.0)
        prony = PronySpectrum(g, lambda)

        @test length(prony.g) == 3
        @test isapprox(sum(prony.g), 1.0, atol=1e-10)
        @test isapprox(relaxation_function(prony, 0.0), 1.0, atol=1e-12)

        val_z1 = 0.2 * exp(-1.0 / 0.1) + 0.5 * exp(-1.0 / 1.0) + 0.3 * exp(-1.0 / 10.0)
        @test isapprox(relaxation_function(prony, 1.0), val_z1, atol=1e-10)

        # Invalid inputs
        @test_throws ArgumentError PronySpectrum(SVector{2, Float64}(0.5, 0.2), SVector{2, Float64}(1.0, 2.0)) # sum != 1
        @test_throws ArgumentError PronySpectrum(SVector{2, Float64}(0.5, 0.5), SVector{2, Float64}(-1.0, 2.0)) # lambda <= 0
    end

    @testset "KWW to Prony Collocation" begin
        kww = KWWSpectrum(0.5)
        prony = kww_to_prony(kww; N=8, z_min=1e-3, z_max=1e3)

        @test isapprox(sum(prony.g), 1.0, atol=1e-6)
        @test all(prony.g .>= 0.0)
        @test all(prony.lambda .> 0.0)

        # Test approximation accuracy across 4 orders of magnitude in z
        test_z = [0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0]
        for z in test_z
            val_kww = relaxation_function(kww, z)
            val_prony = relaxation_function(prony, z)
            @test isapprox(val_kww, val_prony, atol=0.05)
        end
    end
end
