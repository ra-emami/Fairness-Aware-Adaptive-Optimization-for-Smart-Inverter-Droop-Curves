# Fairness-Aware Adaptive Optimization for Smart Inverter Droop Curves

This repository contains code for the adaptive and fairness-aware optimization of smart inverter (SI) Volt–Var (Q–V) droop curves in distribution grids. The optimization model leverages the Current-Voltage AC Optimal Power Flow (IVACOPF) for accurate system modeling and introduces adaptive breakpoints for the Q-V droop curves.

## Key Features:
- **Adaptive Droop Optimization:** The framework updates the breakpoint voltages dynamically based on evolving grid states.
- **Fairness-Aware Control:** Fair allocation schemes (e.g., egalitarian, proportional, priority-based) for PV curtailments are integrated to ensure equitable curtailment distribution across distributed PV units.
- **Enhanced IVACOPF Formulation:** Incorporates power loss constraints to improve model accuracy and computational efficiency.

## Requirements:
- **JuMP**: For optimization modeling.
- **Gurobi** (or alternative solvers like GLPK, Ipopt) for solving the optimization problem.
- **Julia**: Programming language.

## Contributions:

Contributions are welcome! If you'd like to contribute to this project, please feel free to fork the repository, submit issues, and create pull requests. If you have suggestions for improvements or find any bugs, open an issue and we'll work on it together!

Please follow these steps when contributing:

1. Fork the repository.
2. Create a new branch for your feature or bug fix.
3. Ensure the code follows the existing structure and conventions.
4. Submit a pull request describing your changes.

## Contact:

For any questions or feedback, feel free to contact the authors:

- **Rahmat Emami Mirak** - r.emamimirak@gmail.com

## Citation:

If you use this framework in your work, please cite our paper:

> "Fairness-Aware and Adaptive Optimization of Smart Inverter Droop Curves in Distribution Grids," PSCC 2026.
