using Test
using TNM
using Random

@testset "Fitting and Optimization Engines" begin
    # True synthetic parameters
    true_ln_A = -35.0
    true_dh_R = 40000.0
    true_x = 0.5
    true_beta = 0.5
    true_params = [true_ln_A, true_dh_R, true_x, true_beta]

    # Generate synthetic target curve: Cooling 600K -> 450K, Reheating 450K -> 620K at 5 K/s
    proto = MultiRampProtocol(600.0, [450.0, 620.0], [-5.0, 5.0])
    true_model = TNMModel(ArrheniusKinetics(true_ln_A, true_dh_R, true_x), KWWSpectrum(true_beta))
    true_res = solve(true_model, proto, ConvolutionSolver(n_steps=400))
    true_cp = normalized_heat_capacity(true_res)

    exp_data = ExperimentalDSCData(true_res.temperatures, true_cp, proto)
    cost_fn = create_tnm_cost_function([exp_data]; model_type=:arrhenius, solver=ConvolutionSolver(n_steps=300))

    # Test cost at ground truth
    cost_truth = cost_fn(true_params)
    @test cost_truth < 1e-4

    # Optimization parameter bounds
    bounds = [
        (-40.0, -30.0),      # ln_A
        (35000.0, 45000.0),  # dh_R
        (0.3, 0.7),          # x
        (0.3, 0.7)           # beta
    ]

    @testset "Differential Evolution (DE)" begin
        de_res = optimize_differential_evolution(cost_fn, bounds; 
                                                strategy=:best1bin,
                                                max_iters=25, 
                                                pop_multiplier=6, 
                                                top_k=3, 
                                                seed=42)

        @test de_res.best_cost < 0.05
        @test length(de_res.top_k_params) == 3
        @test isapprox(de_res.best_params[3], true_x, atol=0.2)
        @test isapprox(de_res.best_params[4], true_beta, atol=0.2)
    end

    @testset "Particle Swarm Optimization (PSO)" begin
        pso_res = optimize_particle_swarm(cost_fn, bounds; 
                                         max_iters=25, 
                                         n_particles=20, 
                                         top_k=3, 
                                         seed=42)

        @test pso_res.best_cost < 0.08
        @test length(pso_res.top_k_params) == 3
        @test isapprox(pso_res.best_params[3], true_x, atol=0.2)
    end
end
