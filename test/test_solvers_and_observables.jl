using Test
using TNM
using StaticArrays

@testset "Solvers and Observables" begin
    # Reference model parameters (typical oxide / polymer glass)
    ln_A = -75.0
    dh_R = 40000.0
    x = 0.5
    beta = 0.5

    kin = ArrheniusKinetics(ln_A, dh_R, x)
    spec_kww = KWWSpectrum(beta)
    model_kww = TNMModel(kin, spec_kww)

    # Fast cooling and reheating protocol
    # Cool from 600K to 400K at 10 K/s, then reheat from 400K to 650K at 10 K/s
    protocol = MultiRampProtocol(600.0, [400.0, 650.0], [-10.0, 10.0])

    @testset "Convolution Solver (KWW)" begin
        conv_solver = ConvolutionSolver(n_steps=1000)
        res_conv = solve(model_kww, protocol, conv_solver)

        @test length(res_conv.times) == 1000
        @test isapprox(res_conv.fictive_temperatures[1], 600.0, atol=1e-5)
        # Structural freezing: low temperature Tf should freeze above 400K
        @test res_conv.fictive_temperatures[500] > 400.0

        # Observables
        cp_n = normalized_heat_capacity(res_conv)
        @test length(cp_n) == 1000

        # Heating scan peak (overshoot)
        heating_indices = 501:1000
        max_cp = maximum(cp_n[heating_indices])
        @test max_cp > 1.0 # Overshoot characteristic of structural recovery

        # Limiting fictive temperature
        tf_prime = limiting_fictive_temperature(res_conv)
        @test 450.0 < tf_prime < 580.0

        # Glass transition temperature Tg
        tg_peak = glass_transition_temperature(res_conv; method=:peak)
        tg_half = glass_transition_temperature(res_conv; method=:half_step)
        @test 450.0 < tg_peak < 600.0
        @test 450.0 < tg_half < 600.0
    end

    @testset "ODE Solver vs Convolution Solver Agreement" begin
        conv_solver = ConvolutionSolver(n_steps=1500)
        res_conv = solve(model_kww, protocol, conv_solver)

        ode_solver = PronyODESolver(modes=10)
        res_ode = solve(model_kww, protocol, ode_solver)

        @test length(res_ode.times) > 100
        # Check that both solvers predict essentially the same limiting Tf'
        tf_conv = limiting_fictive_temperature(res_conv)
        tf_ode = limiting_fictive_temperature(res_ode)

        @test isapprox(tf_conv, tf_ode, rtol=0.08)
    end

    @testset "TMDSC Analysis" begin
        # TMDSC across glass transition region (480 K to 560 K at 1 K/min, amplitude 0.5 K, period 60 s)
        tmdsc_proto = tmdsc_protocol_k_min(480.0, 560.0, 1.0, 0.5, 60.0)
        conv_solver = ConvolutionSolver(n_steps=2400)
        res_tmdsc = solve(model_kww, tmdsc_proto, conv_solver)

        tmdsc_analysis = analyze_tmdsc(res_tmdsc; Cpg=1.0, delta_Cp=1.0, window_periods=1)
        @test length(tmdsc_analysis.temperatures) > 0
        @test all(tmdsc_analysis.cp_modulus .> 0.0)
        @test all(tmdsc_analysis.cp_real .> 0.0)
        
        # In glassy state (T < 490), |Cp*| ≈ Cpg = 1.0
        # In liquid state (T > 550), |Cp*| ≈ Cpl = 2.0
        @test tmdsc_analysis.cp_modulus[1] < tmdsc_analysis.cp_modulus[end]
    end
end
