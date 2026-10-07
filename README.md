# TNM.jl

[![Julia 1.10+](https://img.shields.io/badge/julia-v1.10+-blue.svg)](https://julialang.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Build Status](https://img.shields.io/badge/tests-119%2F119%20passing-brightgreen.svg)]()

**TNM.jl** is a high-performance, modular, and mathematically rigorous Julia package for modeling, simulating, and evaluating non-linear structural relaxation, physical aging, vitrification, and enthalpy recovery in amorphous materials (glasses, polymers, and supercooled liquids) using the **Tool-Narayanaswamy-Moynihan (TNM)** framework and its modern extensions.

---

## Features

- **Kinetic Formulations**:
  - **Arrhenius (Classical TNM)**: Non-linear activation with apparent activation energy $\Delta h^*/R$ and non-linearity partitioning parameter $x \in (0, 1]$.
  - **VFTH (Vogel-Fulcher-Tammann-Hesse)**: Non-linear relaxation for fragile supercooled liquids.
  - **Adam-Gibbs**: Configurational entropy-based relaxation dynamics.
- **Relaxation Spectra**:
  - **KWW Stretched Exponential**: $\Phi(z) = \exp(-z^\beta)$ with non-exponentiality parameter $0 < \beta \le 1$.
  - **Discrete Prony Series**: $\Phi(z) = \sum_{i=1}^N g_i \exp(-z / \lambda_i)$ with strict normalization $\sum g_i = 1$.
  - **Accelerated NNLS Collocation**: Fast Non-Negative Least Squares collocation using projected gradient descent to discretize continuous KWW spectra into $N$-mode Maxwell representations.
- **Thermal Protocols**:
  - `LinearRampProtocol`: Unidirectional heating and cooling scans.
  - `MultiRampProtocol`: Multi-segment thermal histories.
  - `IsothermalAgingProtocol`: Quench-anneal-reheat cycles for physical aging and structural recovery.
  - `KovacsCrossoverProtocol`: Classical double-temperature jump and isothermal memory effect.
  - `TMDSCProtocol`: Sinusoidally modulated temperature programs $T(t) = T_0 + q_0 t + A_T \sin(\omega t)$.
  - `TabularProtocol`: Ingestion and interpolation of arbitrary experimental time-temperature data.
- **Dual Numerical Solvers**:
  - **Recursive Convolution Solver** (`ConvolutionSolver`): Discrete Moynihan-Markovsky integral algorithm evaluating DSC and aging scans in milliseconds.
  - **State-Space Stiff ODE Solver** (`PronyODESolver`): Multi-mode differential equation formulation solved with high-order stiff algorithms from `OrdinaryDiffEq.jl` (`RadauIIA5`, `Rodas5P`).
- **Physical Observables & Analysis**:
  - Limiting fictive temperature $T_f'$ via Moynihan equal-area construction.
  - Glass transition temperature $T_g$ determination (`:peak`, `:half_step`, `:inflection`).
  - Apparent and normalized heat capacity ($C_p(T)$, $C_p^N = dT_f/dT$).
  - Sliding Fourier transform deconvolution for TMDSC (storage capacity $C'$, loss capacity $C''$, complex modulus $|C_p^*|$, and phase angle $\delta$).
- **Metaheuristic Parameter Optimization**:
  - **Differential Evolution** (`optimize_differential_evolution`): `DE/best1bin` and `DE/rand1bin` with mutation dither, binomial crossover, parameter bound clamping, and top-$K$ candidate archiving.
  - **Particle Swarm Optimization** (`optimize_particle_swarm`): Multi-particle swarm with inertia decay, cognitive/social velocity acceleration, velocity clamping, and top-$K$ tracking.
  - Simultaneous multi-rate / multi-protocol fitting against experimental DSC data.
- **Unit Systems**:
  - Standard physical units (`PhysicalUnits` in Kelvin and seconds).
  - Dimensionless reduced units (`ReducedUnits`) for master-curve analysis ($\theta = T/T_{\text{ref}}$, $\tilde{t} = t/\tau_{\text{ref}}$, $\tilde{q} = q\tau_{\text{ref}}/T_{\text{ref}}$, $\epsilon = \Delta h^*/(RT_{\text{ref}})$).
- **Data Export & Gnuplot Integration**:
  - Clean whitespace-delimited `.dat` files with `# prop1 prop2 ...` column headers.
  - Ready-to-run `.gnu` script generation for instantaneous plotting with Gnuplot.

---

## Installation

Clone the repository and instantiate the environment:

```bash
git clone https://github.com/Riperedo/TNM.jl.git
cd TNM.jl
julia --project=. -e 'import Pkg; Pkg.instantiate()'
```

---

## Quick Start

### 1. Classical DSC Heating and Cooling Scan

```julia
using TNM

# 1. Define Kinetics and Spectrum
kinetics = ArrheniusKinetics(-75.0, 42000.0, 0.5)  # ln(A), Δh*/R, x
spectrum = KWWSpectrum(0.6)                         # β
model = TNMModel(kinetics, spectrum)

# 2. Thermal Protocol: Cool 650 K -> 450 K at 10 K/min, Reheat 450 K -> 650 K at 10 K/min
q = 10.0 / 60.0  # K/s
protocol = MultiRampProtocol(650.0, [450.0, 650.0], [-q, q])

# 3. Solve
solver = ConvolutionSolver(n_steps=2000)
result = solve(model, protocol, solver)

# 4. Extract Physical Observables
tf_prime = limiting_fictive_temperature(result)
tg_peak  = glass_transition_temperature(result; method=:peak)
println("Limiting Tf': ", tf_prime, " K")
println("Glass Transition Tg: ", tg_peak, " K")

# 5. Export Data and Gnuplot Script
export_dat("output/01_classical_tnm_dsc.dat", result; Cpg=1.0, Cpl=2.0)
generate_gnuplot_script("output/01_classical_tnm_dsc.dat", "output/01_classical_tnm_dsc.gnu";
                        title="TNM DSC Scan", plot_type=:cp_normalized)
```

To plot the exported data directly:
```bash
gnuplot -p output/01_classical_tnm_dsc.gnu
```

---

### 2. Kovacs Crossover Memory Effect

```julia
using TNM

kinetics = ArrheniusKinetics(-75.0, 40000.0, 0.4)
spectrum = kww_to_prony(KWWSpectrum(0.5); N=8)
model = TNMModel(kinetics, spectrum)

# Quench from 550 K to 470 K (hold 300 s), jump to 490 K (observe 2000 s)
protocol = KovacsCrossoverProtocol(550.0, 470.0, 490.0, 300.0, 2000.0; q_quench=50.0)

result = solve(model, protocol, ConvolutionSolver(n_steps=3000))
export_dat("output/02_kovacs_crossover.dat", result)
generate_gnuplot_script("output/02_kovacs_crossover.dat", "output/02_kovacs_crossover.gnu";
                        title="Kovacs Memory Effect", plot_type=:fictive_temperature)
```

---

### 3. Temperature-Modulated DSC (TMDSC)

```julia
using TNM

model = TNMModel(ArrheniusKinetics(-75.0, 42000.0, 0.5), KWWSpectrum(0.55))

# Modulated ramp: 480 K to 560 K at 1 K/min, Amplitude = 0.5 K, Period = 60 s
protocol = tmdsc_protocol_k_min(480.0, 560.0, 1.0, 0.5, 60.0)

# Simulate and perform sliding Fourier deconvolution
res = solve(model, protocol, ConvolutionSolver(n_steps=2400))
tmdsc_res = analyze_tmdsc(res; Cpg=1.0, delta_Cp=1.0, window_periods=1)

export_dat("output/03_tmdsc_analysis.dat", tmdsc_res)
generate_gnuplot_script("output/03_tmdsc_analysis.dat", "output/03_tmdsc_analysis.gnu";
                        title="TMDSC Analysis", plot_type=:tmdsc)
```

---

### 4. Parameter Estimation (DE & PSO)

```julia
using TNM

# Ingest multi-rate experimental curves
exp1 = load_experimental_dsc("data/scan_10k_min.dat", protocol1)
exp2 = load_experimental_dsc("data/scan_20k_min.dat", protocol2)

# Build objective function
cost_fn = create_tnm_cost_function([exp1, exp2]; model_type=:arrhenius)

# Parameter search bounds: [ln_A, Δh*/R, x, β]
bounds = [
    (-85.0, -65.0),      # ln_A
    (35000.0, 50000.0),  # Δh*/R
    (0.2, 0.8),          # x
    (0.3, 0.8)           # β
]

# Run Differential Evolution with Top-K Archiving
de_res = optimize_differential_evolution(cost_fn, bounds;
                                        strategy=:best1bin,
                                        max_iters=35,
                                        top_k=5,
                                        verbose=true)

println("Best Cost: ", de_res.best_cost)
println("Recovered Parameters: ", de_res.best_params)
println("Top 5 Candidate Solutions: ", de_res.top_k_params)
```

---

## Running Examples and Tests

Run the full automated test suite (119 unit and integration tests):
```bash
julia --project=. test/runtests.jl
```

Run individual example scripts:
```bash
julia --project=. examples/01_classical_tnm_dsc.jl
julia --project=. examples/02_kovacs_memory_effect.jl
julia --project=. examples/03_weyer_tmdsc.jl
julia --project=. examples/04_parameter_fitting.jl
```

---

## Repository Structure

```
TNM.jl/
├── Project.toml
├── README.md
├── src/
│   ├── TNM.jl                         # Main module and public exports
│   ├── types.jl                       # Abstract types and unit systems
│   ├── kinetics/                      # Relaxation kinetics (Arrhenius, VFTH, Adam-Gibbs)
│   ├── spectrum/                      # KWW, Prony, and NNLS collocation
│   ├── protocols/                     # Thermal protocols (Ramps, Aging, Kovacs, TMDSC, Tabular)
│   ├── solvers/                       # Moynihan convolution & stiff ODE solvers
│   ├── observables/                   # Fictive temperature Tf', Tg, CpN, TMDSC analysis
│   ├── export/                        # .dat formatter & Gnuplot script generator
│   ├── fitting/                       # DE and PSO metaheuristic optimization engines
│   └── data/                          # Experimental data ingestion and loaders
├── examples/                          # 4 complete runnable workflow examples
└── test/                              # Comprehensive test suites (119 tests)
```

---

## License

This project is licensed under the [MIT License](LICENSE).
