# Example 03: Temperature-Modulated DSC (TMDSC) Simulation
# Simulates sinusoidally modulated DSC across the glass transition and extracts complex heat capacity C*, C', C'', and phase lag.

using TNM
using Printf

function run_example()
    println("="^70)
    println("TNM.jl - Example 03: Temperature-Modulated DSC (TMDSC)")
    println("="^70)

    # 1. Model setup
    ln_A = -75.0
    dh_R = 42000.0
    x = 0.5
    beta = 0.55
    model = TNMModel(ArrheniusKinetics(ln_A, dh_R, x), KWWSpectrum(beta))

    # 2. TMDSC Protocol: Underlying ramp 480 K -> 560 K at 1 K/min with modulation:
    # Amplitude A_T = 0.5 K, Period t_p = 60 s
    protocol = tmdsc_protocol_k_min(480.0, 560.0, 1.0, 0.5, 60.0)

    println("TMDSC Protocol Parameters:")
    println("  Temperature Range: 480 K to 560 K")
    println("  Underlying Rate:   1 K/min")
    println("  Modulation Amp:    0.5 K")
    println("  Modulation Period: 60 s")
    println("  Total Duration:    $(duration(protocol)) s")

    # 3. Simulate high-resolution time series
    println("\nRunning high-resolution TMDSC simulation...")
    # ~30 points per 60s period -> ~2400 steps
    solver = ConvolutionSolver(n_steps=2400)
    res = solve(model, protocol, solver)

    # 4. Perform sliding Fourier deconvolution
    println("Performing sliding Fourier transform analysis...")
    tmdsc_res = analyze_tmdsc(res; Cpg=1.0, delta_Cp=1.0, window_periods=1)

    println("\nAnalysis Results Summary:")
    println("  Evaluated sliding points: $(length(tmdsc_res.times))")
    @printf("  Glassy Modulus |Cp*|:   %.4f (T = %.1f K)\n", tmdsc_res.cp_modulus[1], tmdsc_res.temperatures[1])
    @printf("  Liquid Modulus |Cp*|:   %.4f (T = %.1f K)\n", tmdsc_res.cp_modulus[end], tmdsc_res.temperatures[end])
    @printf("  Peak Phase Lag δ_max:   %.2f deg\n", maximum(tmdsc_res.phase_deg))

    # 5. Export results
    output_dir = joinpath(@__DIR__, "..", "output")
    mkpath(output_dir)
    dat_file = joinpath(output_dir, "03_tmdsc_analysis.dat")
    gnu_file = joinpath(output_dir, "03_tmdsc_analysis.gnu")

    export_dat(dat_file, tmdsc_res)
    generate_gnuplot_script(dat_file, gnu_file; title="TMDSC Complex Heat Capacity and Phase", plot_type=:tmdsc)

    println("\nGenerated Output Files:")
    println("  Data file:     $(dat_file)")
    println("  Gnuplot file:  $(gnu_file)")
    println("  To plot: run `gnuplot -p $(gnu_file)`")
    println("="^70)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_example()
end
