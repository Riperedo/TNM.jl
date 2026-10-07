# Example 02: Kovacs Crossover / Memory Effect Simulation
# Demonstrates non-monotonic structural relaxation resulting from the multi-mode relaxation spectrum.

using TNM
using Printf

function run_example()
    println("="^70)
    println("TNM.jl - Example 02: Kovacs Crossover Memory Effect")
    println("="^70)

    # 1. Multi-mode Prony series spectrum (captures broad distribution of relaxation times)
    ln_A = -75.0
    dh_R = 40000.0
    x = 0.4
    kinetics = ArrheniusKinetics(ln_A, dh_R, x)
    
    # 8-mode Prony representation of stretched exponential beta = 0.5
    spectrum = kww_to_prony(KWWSpectrum(0.5); N=8)
    model = TNMModel(kinetics, spectrum)

    # 2. Kovacs Crossover Protocol:
    # 1) Quench from equilibrium T0 = 550 K to T_low = 470 K
    # 2) Anneal at 470 K for t_anneal = 300 s
    # 3) Jump to T_target = 490 K and observe isothermal relaxation for 2000 s
    protocol = KovacsCrossoverProtocol(550.0, 470.0, 490.0, 300.0, 2000.0; q_quench=50.0)

    println("Protocol parameters:")
    println("  Initial Equilibrium T0:  550 K")
    println("  Annealing Temp T_low:    470 K (300 s)")
    println("  Observation Temp T_target: 490 K (2000 s)")
    println("  Total Duration:          $(duration(protocol)) s")

    # 3. Solve using Prony Convolution Solver
    println("\nSolving structural relaxation trajectory...")
    solver = ConvolutionSolver(n_steps=3000)
    res = solve(model, protocol, solver)

    # 4. Export
    output_dir = joinpath(@__DIR__, "..", "output")
    mkpath(output_dir)
    dat_file = joinpath(output_dir, "02_kovacs_crossover.dat")
    gnu_file = joinpath(output_dir, "02_kovacs_crossover.gnu")

    export_dat(dat_file, res)
    generate_gnuplot_script(dat_file, gnu_file; title="Kovacs Crossover Memory Effect", plot_type=:fictive_temperature)

    println("\nGenerated Output Files:")
    println("  Data file:     $(dat_file)")
    println("  Gnuplot file:  $(gnu_file)")
    println("  To plot: run `gnuplot -p $(gnu_file)`")
    println("="^70)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_example()
end
