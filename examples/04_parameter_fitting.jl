# Example 04: Parameter Estimation via Metaheuristics (DE and PSO)
# Recovers TNM parameters from multi-rate synthetic experimental DSC curves.

using TNM
using Printf

function run_example()
    println("="^70)
    println("TNM.jl - Example 04: Parameter Estimation with DE and PSO")
    println("="^70)

    # 1. Ground Truth Parameters
    true_ln_A = -75.0
    true_dh_R = 42000.0
    true_x = 0.52
    true_beta = 0.58
    true_params = [true_ln_A, true_dh_R, true_x, true_beta]

    println("Ground Truth Model Parameters:")
    println("  ln(A)  = $(true_ln_A)")
    println("  Δh*/R  = $(true_dh_R) K")
    println("  x      = $(true_x)")
    println("  β      = $(true_beta)")

    true_kin = ArrheniusKinetics(true_ln_A, true_dh_R, true_x)
    true_model = TNMModel(true_kin, KWWSpectrum(true_beta))

    # 2. Generate Synthetic Multi-Rate Experimental Curves (Cool at 20 K/min, reheat at 10 K/min & 20 K/min)
    println("\nGenerating synthetic multi-rate experimental curves...")
    proto1 = MultiRampProtocol(620.0, [460.0, 620.0], [-20.0 / 60.0, 10.0 / 60.0])
    proto2 = MultiRampProtocol(620.0, [460.0, 620.0], [-20.0 / 60.0, 20.0 / 60.0])

    res1 = solve(true_model, proto1, ConvolutionSolver(n_steps=400))
    res2 = solve(true_model, proto2, ConvolutionSolver(n_steps=400))

    exp1 = ExperimentalDSCData(res1.temperatures, normalized_heat_capacity(res1), proto1)
    exp2 = ExperimentalDSCData(res2.temperatures, normalized_heat_capacity(res2), proto2)

    # 3. Formulate Combined Objective Function
    cost_fn = create_tnm_cost_function([exp1, exp2]; model_type=:arrhenius, solver=ConvolutionSolver(n_steps=300))

    # Search bounds: [ln_A, dh_R, x, beta]
    bounds = [
        (-85.0, -65.0),      # ln_A
        (35000.0, 50000.0),  # dh_R
        (0.2, 0.8),          # x
        (0.3, 0.8)           # beta
    ]

    # 4. Optimize using Differential Evolution (best1bin)
    println("\n--- Running Differential Evolution (best1bin) ---")
    de_res = optimize_differential_evolution(cost_fn, bounds;
                                            strategy=:best1bin,
                                            max_iters=35,
                                            pop_multiplier=6,
                                            top_k=3,
                                            seed=42,
                                            verbose=true)

    println("\nDE Optimization Complete:")
    println("  Best Cost (MSE): $(de_res.best_cost)")
    println("  Recovered Parameters: [ln_A, Δh*/R, x, β]")
    @printf("    ln(A) = %.4f (true: %.4f)\n", de_res.best_params[1], true_ln_A)
    @printf("    Δh*/R = %.1f (true: %.1f)\n", de_res.best_params[2], true_dh_R)
    @printf("    x     = %.4f (true: %.4f)\n", de_res.best_params[3], true_x)
    @printf("    β     = %.4f (true: %.4f)\n", de_res.best_params[4], true_beta)

    # 5. Optimize using Particle Swarm Optimization (PSO)
    println("\n--- Running Particle Swarm Optimization (PSO) ---")
    pso_res = optimize_particle_swarm(cost_fn, bounds;
                                      max_iters=35,
                                      n_particles=24,
                                      top_k=3,
                                      seed=42,
                                      verbose=true)

    println("\nPSO Optimization Complete:")
    println("  Best Cost (MSE): $(pso_res.best_cost)")
    println("  Recovered Parameters: [ln_A, Δh*/R, x, β]")
    @printf("    ln(A) = %.4f (true: %.4f)\n", pso_res.best_params[1], true_ln_A)
    @printf("    Δh*/R = %.1f (true: %.1f)\n", pso_res.best_params[2], true_dh_R)
    @printf("    x     = %.4f (true: %.4f)\n", pso_res.best_params[3], true_x)
    @printf("    β     = %.4f (true: %.4f)\n", pso_res.best_params[4], true_beta)

    println("="^70)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_example()
end
