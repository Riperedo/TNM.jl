using Test
using TNM

@testset "Units and Kinetics" begin
    @testset "Unit Systems" begin
        phys = PhysicalUnits()
        @test phys isa AbstractUnitSystem

        red = ReducedUnits(600.0, 10.0)
        @test red isa AbstractUnitSystem
        @test red.T_ref == 600.0
        @test red.tau_ref == 10.0

        @test_throws ArgumentError ReducedUnits(-100.0)
        @test_throws ArgumentError ReducedUnits(600.0, -1.0)
    end

    @testset "Arrhenius Kinetics" begin
        ln_A = -35.0
        dh_R = 50000.0
        x = 0.5
        kin = ArrheniusKinetics(ln_A, dh_R, x)

        @test kin.ln_A == ln_A
        @test kin.dh_R == dh_R
        @test kin.x == x

        # Temperature = Tf = 600 K (equilibrium)
        T_eq = 600.0
        ln_tau = log_relaxation_time(kin, T_eq, T_eq, PhysicalUnits())
        expected_ln_tau = -35.0 + 50000.0 * (0.5 / 600.0 + 0.5 / 600.0)
        @test isapprox(ln_tau, expected_ln_tau, atol=1e-10)

        tau = relaxation_time(kin, T_eq, T_eq, PhysicalUnits())
        @test isapprox(tau, exp(expected_ln_tau), rtol=1e-10)

        # Non-equilibrium T = 550, Tf = 600
        ln_tau_neq = log_relaxation_time(kin, 550.0, 600.0, PhysicalUnits())
        expected_neq = -35.0 + 50000.0 * (0.5 / 550.0 + 0.5 / 600.0)
        @test isapprox(ln_tau_neq, expected_neq, atol=1e-10)

        # Reduced units test: θ = 1.0, θ_f = 1.0 at T_ref = 600 K, tau_ref = exp(expected_ln_tau)
        tau_ref_val = exp(expected_ln_tau)
        red_units = ReducedUnits(600.0, tau_ref_val)
        ln_tau_red = log_relaxation_time(kin, 1.0, 1.0, red_units)
        @test isapprox(ln_tau_red, 0.0, atol=1e-10)

        # Invalid bounds
        @test_throws ArgumentError ArrheniusKinetics(ln_A, -100.0, x)
        @test_throws ArgumentError ArrheniusKinetics(ln_A, dh_R, 0.0)
        @test_throws ArgumentError ArrheniusKinetics(ln_A, dh_R, 1.2)
        @test_throws ArgumentError log_relaxation_time(kin, -300.0, 600.0, PhysicalUnits())
    end

    @testset "VFTH Kinetics" begin
        ln_A = -30.0
        B = 3000.0
        T_inf = 350.0
        x = 0.6
        vfth = VFTHKinetics(ln_A, B, T_inf, x)

        @test vfth.B == B
        @test vfth.T_inf == T_inf

        T = 400.0
        Tf = 410.0
        ln_tau = log_relaxation_time(vfth, T, Tf, PhysicalUnits())
        expected = -30.0 + 3000.0 * (0.6 / (400.0 - 350.0) + 0.4 / (410.0 - 350.0))
        @test isapprox(ln_tau, expected, atol=1e-10)

        # Divergence / error below T_inf
        @test_throws ArgumentError log_relaxation_time(vfth, 340.0, 410.0, PhysicalUnits())
        @test_throws ArgumentError log_relaxation_time(vfth, 400.0, 340.0, PhysicalUnits())
    end

    @testset "Adams-Gibbs Kinetics" begin
        ln_A = -30.0
        C = 5000.0
        Sc0 = 10.0
        delta_Cp = 0.5
        T0 = 400.0
        ag = AdamsGibbsKinetics(ln_A, C, Sc0, delta_Cp, T0)

        @test ag.C == C
        @test ag.Sc0 == Sc0

        T = 400.0
        Tf = 400.0
        ln_tau = log_relaxation_time(ag, T, Tf, PhysicalUnits())
        expected = -30.0 + 5000.0 / (400.0 * 10.0)
        @test isapprox(ln_tau, expected, atol=1e-10)

        @test_throws ArgumentError AdamsGibbsKinetics(ln_A, -10.0, Sc0, delta_Cp, T0)
    end
end
