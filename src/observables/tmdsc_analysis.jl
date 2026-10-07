# Observables: Temperature-Modulated DSC (TMDSC) Complex Heat Capacity Analysis

using Statistics

"""
    TMDSCResult

Container for Temperature-Modulated DSC analysis.

# Fields
- `times::Vector{Float64}`: Center times for sliding windows.
- `temperatures::Vector{Float64}`: Center underlying temperatures.
- `cp_real::Vector{Float64}`: Real part of complex heat capacity C' (in-phase / storage).
- `cp_imag::Vector{Float64}`: Imaginary part of complex heat capacity C'' (out-of-phase / loss).
- `cp_modulus::Vector{Float64}`: Modulus of complex heat capacity |Cp*| (reversing Cp).
- `phase_deg::Vector{Float64}`: Phase angle δ_S in degrees.
"""
struct TMDSCResult
    times::Vector{Float64}
    temperatures::Vector{Float64}
    cp_real::Vector{Float64}
    cp_imag::Vector{Float64}
    cp_modulus::Vector{Float64}
    phase_deg::Vector{Float64}
end

"""
    analyze_tmdsc(res::SimulationResult; Cpg::Real=1.0, delta_Cp::Real=1.0, window_periods::Int=1)

Extracts the complex heat capacity components (C', C'', |Cp*|, and phase angle δ) across sliding modulation periods
from a TMDSC simulation result.
"""
function analyze_tmdsc(res::SimulationResult; Cpg::Real=1.0, delta_Cp::Real=1.0, window_periods::Int=1)
    res.protocol isa TMDSCProtocol || throw(ArgumentError("Protocol must be a TMDSCProtocol"))
    
    proto = res.protocol
    period = proto.period
    omega = proto.omega
    t = res.times
    T = res.temperatures
    hf = heat_flow(res; Cpg=Cpg, delta_Cp=delta_Cp)
    q = res.heating_rates

    dt_step = t[2] - t[1]
    pts_per_period = round(Int, period / dt_step)
    window_pts = pts_per_period * window_periods
    
    N = length(t)
    N > window_pts || throw(ArgumentError("Not enough time points for TMDSC analysis"))

    stride = max(1, div(pts_per_period, 4))
    
    out_times = Float64[]
    out_T = Float64[]
    out_c_real = Float64[]
    out_c_imag = Float64[]
    out_modulus = Float64[]
    out_phase = Float64[]

    for start_idx in 1:stride:(N - window_pts)
        end_idx = start_idx + window_pts - 1
        t_win = @view t[start_idx:end_idx]
        T_win = @view T[start_idx:end_idx]
        q_win = @view q[start_idx:end_idx]
        hf_win = @view hf[start_idx:end_idx]

        # Center time and temperature
        center_t = mean(t_win)
        center_T = mean(T_win)

        # Discrete Fourier Transform at fundamental frequency ω
        # A* = (2/T_p) * ∫ f(t) * exp(-i ω t) dt
        cos_wt = cos.(omega .* t_win)
        sin_wt = sin.(omega .* t_win)

        # Heating rate complex amplitude
        Aq_cos = mean(q_win .* cos_wt) * 2.0
        Aq_sin = -mean(q_win .* sin_wt) * 2.0
        Aq_complex = Complex(Aq_cos, Aq_sin)

        # Heat flow complex amplitude
        Ahf_cos = mean(hf_win .* cos_wt) * 2.0
        Ahf_sin = -mean(hf_win .* sin_wt) * 2.0
        Ahf_complex = Complex(Ahf_cos, Ahf_sin)

        if abs(Aq_complex) > 1e-12
            Cp_complex = Ahf_complex / Aq_complex
            c_real = real(Cp_complex)
            c_imag = -imag(Cp_complex) # Standard convention in thermal physics: Cp* = Cp' - i Cp''
            modulus = abs(Cp_complex)
            phase = rad2deg(atan(c_imag, c_real))

            push!(out_times, center_t)
            push!(out_T, center_T)
            push!(out_c_real, c_real)
            push!(out_c_imag, c_imag)
            push!(out_modulus, modulus)
            push!(out_phase, phase)
        end
    end

    return TMDSCResult(out_times, out_T, out_c_real, out_c_imag, out_modulus, out_phase)
end
