# Example 01: Classical Tool-Narayanaswamy-Moynihan DSC Simulation
# Simulates standard cooling followed by reheating across the glass transition.

using TNM
using Printf

function run_example()
    println("="^70)
    println("TNM.jl - Example 01: Classical Cooling and Heating Scan")
    println("="^70)

    # 1. Define Kinetics and Spectrum (Representative parameters for oxide/polymeric glass)
    ln_A = -75.0          # Pre-exponential factor ln(A / s)
    dh_R = 42000.0        # Apparent activation energy Δh*/R (K)
    x = 0.5               # Non-linearity parameter
    beta = 0.6            # KWW stretched exponential exponent

    kinetics = ArrheniusKinetics(ln_A, dh_R, x)
    spectrum = KWWSpectrum(beta)
    model = TNMModel(kinetics, spectrum)

    println("Model Parameters:")
    println("  ln(A)  = $(ln_A)")
    println("  Δh*/R  = $(dh_R) K")
    println("  x      = $(x)")
    println("  β      = $(beta)")

    # 2. Define Thermal Protocol: Cool from 650 K to 450 K at 10 K/min, then reheat to 650 K at 10 K/min
    q_rate = 10.0 / 60.0  # 10 K/min in K/s
    protocol = MultiRampProtocol(650.0, [450.0, 650.0], [-q_rate, q_rate])
    
    println("\nProtocol: Cool 650 K -> 450 K at 10 K/min, Reheat 450 K -> 650 K at 10 K/min")
    println("Total duration: $(duration(protocol)) s ($(duration(protocol)/60) min)")

    # 3. Solve using Recursive Convolution Solver
    println("\nRunning ConvolutionSolver...")
    solver = ConvolutionSolver(n_steps=2000)
    res = solve(model, protocol, solver)

    # 4. Extract Physical Observables
    tf_prime = limiting_fictive_temperature(res)
    tg_peak = glass_transition_temperature(res; method=:peak)
    tg_half = glass_transition_temperature(res; method=:half_step)

    @printf("\nCalculated Observables:\n")
    @printf("  Limiting Fictive Temperature Tf': %.2f K\n", tf_prime)
    @printf("  Glass Transition Temp Tg (Peak):  %.2f K\n", tg_peak)
    @printf("  Glass Transition Temp Tg (1/2 Cp): %.2f K\n", tg_half)

    # 5. Export to .dat file with # prop1 prop2 ... header
    output_dir = joinpath(@__DIR__, "..", "output")
    mkpath(output_dir)
    dat_file = joinpath(output_dir, "01_classical_tnm_dsc.dat")
    gnu_file = joinpath(output_dir, "01_classical_tnm_dsc.gnu")

    export_dat(dat_file, res; Cpg=1.0, Cpl=2.0, delta_Cp=1.0)
    generate_gnuplot_script(dat_file, gnu_file; title="TNM DSC Scan: Cool & Heat at 10 K/min", plot_type=:cp_normalized)

    println("\nGenerated Output Files:")
    println("  Data file:     $(dat_file)")
    println("  Gnuplot file:  $(gnu_file)")
    println("  To plot: run `gnuplot -p $(gnu_file)`")
    println("="^70)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_example()
end
