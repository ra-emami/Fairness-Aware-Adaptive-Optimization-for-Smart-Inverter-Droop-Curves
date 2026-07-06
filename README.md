# Fairness-Aware Adaptive Optimization of Smart Inverter Droop Curves

[![DOI](https://img.shields.io/badge/DOI-10.1016%2Fj.epsr.2026.113613-blue)](https://doi.org/10.1016/j.epsr.2026.113613)
[![Julia](https://img.shields.io/badge/Julia-1.6%2B-9558B2?logo=julia)](https://julialang.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

Reference implementation for the paper:

> R. Emami Mirak and A. Inaolaji, **"Adaptive and fair optimization of smart inverter droop curves in distribution grids,"** *Electric Power Systems Research*, vol. 262, 2027, Art. no. 113613. [https://doi.org/10.1016/j.epsr.2026.113613](https://doi.org/10.1016/j.epsr.2026.113613)
>
> Presented at the 24th Power Systems Computation Conference (PSCC 2026), Limassol, Cyprus.

The framework determines standards-compliant smart inverter (SI) Volt–Var (Q–V) droop settings with **variable breakpoint voltages** inside a multi-period distribution optimal power flow (DOPF). It combines:

- **Enhanced IVACOPF formulation** — the Current–Voltage AC OPF is augmented with sending-end power-balance and linearized branch-loss constraints, improving numerical accuracy and convergence (maximum linearization error below 10⁻⁶ p.u. within 3–4 iterations).
- **Adaptive droop optimization** — the interior breakpoint voltages of each inverter's Q–V curve are decision variables, re-optimized hourly within IEEE 1547 bounds while load and PV profiles vary at 15-minute resolution.
- **Fairness-aware curtailment** — egalitarian, proportional, and priority-based PV curtailment policies are embedded in the objective; the priority-based scheme is reformulated with big-M constraints so it is directly usable by off-the-shelf MILP/MINLP solvers.

## Repository structure

```
├── data/
│   ├── 33_bus_system_data.csv     # 33-bus feeder branch and peak-load data (12.66 kV, 10 MVA base)
│   ├── load_variation_15min.csv   # 15-min P/Q profiles (%) per load class + bus-to-class mapping
│   └── solar_profile_15min.csv    # 15-min PV generation profile (%)
├── scenarios/
│   ├── scenario1_fixed_bp_no_fairness.jl      # S1: fixed breakpoints, max PV utilization
│   ├── scenario2_fixed_bp_proportional.jl     # S2: fixed breakpoints, proportional fairness
│   ├── scenario3_variable_bp_no_fairness.jl   # S3: variable breakpoints, max PV utilization
│   ├── scenario4_variable_bp_egalitarian.jl   # S4: variable breakpoints, egalitarian fairness
│   ├── scenario5_variable_bp_proportional.jl  # S5: variable breakpoints, proportional fairness
│   └── scenario6_variable_bp_priority.jl      # S6: variable breakpoints, priority-based fairness
├── pareto/
│   ├── pareto_scenario4_egalitarian.jl        # Pareto point for S4 (set α, re-run)
│   ├── pareto_scenario5_proportional.jl       # Pareto point for S5 (set α, re-run)
│   └── pareto_scenario6_priority.jl           # Pareto point for S6 (Fig. 9 of the paper)
├── Project.toml
├── CITATION.cff
├── LICENSE
└── README.md
```

## Case-study scenarios

| Script | Breakpoints | Fairness policy | Paper Eq. |
|---|---|---|---|
| `scenario1_fixed_bp_no_fairness.jl` | fixed {0.88, 0.90, 0.97, 1.00, 1.02, 1.10} | none (min total curtailment) | (36) |
| `scenario2_fixed_bp_proportional.jl` | fixed {0.88, 0.90, 0.97, 1.00, 1.02, 1.10} | proportional | (38) |
| `scenario3_variable_bp_no_fairness.jl` | variable, hourly (Eqs. 25–27) | none (min total curtailment) | (36) |
| `scenario4_variable_bp_egalitarian.jl` | variable, hourly (Eqs. 25–27) | egalitarian | (37) |
| `scenario5_variable_bp_proportional.jl` | variable, hourly (Eqs. 25–27) | proportional | (38) |
| `scenario6_variable_bp_priority.jl` | variable, hourly (Eqs. 25–27) | priority-based (big-M) | (37), (41)–(45) |

All scenarios use the objective weights α = β = 0.5 (Eq. 33). In this implementation α weights the normalized fairness/curtailment term Ψ and β the normalized voltage-deviation term, matching the (α, β) labeling of the Pareto front in Fig. 9 of the paper. The normalization constants η (named `f1`, `f2` in the scripts) are the respective single-objective optima.

## Model summary

- **Grid model:** single-phase Current–Voltage ACOPF (IVACOPF) with first-order Taylor linearization of the voltage-magnitude and power-balance equations around the previous iterate, solved iteratively (Eqs. 1–8). The *enhanced* formulation adds sending-end power balances (Eqs. 9–12) and linearized branch-loss constraints (Eqs. 15–16).
- **Convergence:** the maximum linearization error (MLE) — the largest absolute mismatch between the bilinear branch loss and its linearization — must fall below 10⁻⁶ p.u.; all scenarios converge within 3–4 iterations.
- **SI capability:** 32-sided polygonal approximation of the apparent-power semicircle (Eq. 17) plus active/reactive output limits (Eq. 18).
- **Volt–Var droop:** 5-segment piecewise-linear Q–V curve with six breakpoints and fixed ordinates (+Qmax, +Qmax, 0, 0, −Qmax, −Qmax); the Lambda interpolation method with an SOS2 structure enforced by binary adjacency indicators (Eqs. 19–24). With variable breakpoints the interpolation is bilinear, making the overall problem a MINLP handled by Gurobi's nonconvex machinery.
- **Adaptive breakpoints:** extreme points fixed at 0.80 / 1.20 p.u. (Eq. 27); interior points bounded per IEEE 1547 (Eq. 26): 0.82 ≤ V₂ ≤ V₃ − 0.02, 0.97 ≤ V₃ ≤ 1.00, 1.00 ≤ V₄ ≤ 1.03, V₄ + 0.02 ≤ V₅ ≤ 1.18.

## Test system

Modified 33-bus distribution feeder (12.66 kV, 10 MVA base) with three PV units at buses **7, 18, and 33**, rated **2.2, 1.0, and 1.5 MW** with inverter apparent-power ratings oversized by 10%. Loads are grouped into industrial (buses 26–33), commercial (19–25), and residential (2–18) classes, each following a representative 15-minute active/reactive profile over 24 hours; PV availability follows a 15-minute clear-sky irradiance profile.

Data notes:

- `33_bus_system_data.csv` — row 2 carries the base values (columns 8–9: substation voltage in kV, base power in kVA) alongside the branch table; the scripts skip it accordingly.
- `load_variation_15min.csv` — 96 data rows of six percentage columns (P and Q for the three classes); columns 10–11 hold the bus-to-class mapping; two trailing summary rows ("Average", "Max") are not read by the scripts.
- `solar_profile_15min.csv` — 96 values in percent of rated PV output.

## Requirements

- [Julia](https://julialang.org/) ≥ 1.6
- [Gurobi](https://www.gurobi.com/) with a valid license (free academic licenses are available) and the [Gurobi.jl](https://github.com/jump-dev/Gurobi.jl) wrapper
- Julia packages: JuMP, CSV, DataFrames, Plots (declared in `Project.toml`)

Install the dependencies with:

```bash
julia --project=. -e "using Pkg; Pkg.instantiate()"
```

## Running

Run any scenario from the repository root:

```bash
julia --project=. scenarios/scenario5_variable_bp_proportional.jl
```

Each script iterates the enhanced IVACOPF to convergence and prints, per iteration: solver size and timing, substation/PV energy balances, PV curtailment (total, per unit, and the PVC percentage of Eq. 35), the fairness term, the average voltage deviation (AVD, Eq. 34), and the MLE. Two plots are displayed at the end (PV dispatch/curtailment over 24 h and the network voltage-range profile).

**Expected runtime:** the fixed-breakpoint scenarios (S1–S2) solve in seconds; the variable-breakpoint scenarios (S3–S6) are bilinear MINLPs and complete in under four minutes on a modest laptop (12th Gen Intel i3, 12 GB RAM), consistent with the timings reported in Table 1 of the paper.

**Pareto fronts:** the scripts in `pareto/` compute one trade-off point per run. Set the weight `α` in the script to each value in {0.1, …, 0.9} and re-run to trace the front; Fig. 9 of the paper shows the resulting front for Scenario 6.

## Citation

If you use this code or data in your work, please cite:

```bibtex
@article{EmamiMirak2027adaptive,
  title   = {Adaptive and fair optimization of smart inverter droop curves in distribution grids},
  author  = {Emami Mirak, Rahmat and Inaolaji, Adedoyin},
  journal = {Electric Power Systems Research},
  volume  = {262},
  pages   = {113613},
  year    = {2027},
  doi     = {10.1016/j.epsr.2026.113613}
}
```

## Contributing

Contributions are welcome — feel free to open an issue or submit a pull request.

## Contact

**Rahmat Emami Mirak** — r.emamimirak@gmail.com · ra.emami@aut.ac.ir

## License

Released under the [MIT License](LICENSE).
