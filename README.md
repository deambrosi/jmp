# Dynamic Migration Network Model

This repository contains MATLAB code for simulating a dynamic migration model with network-based help. The code was developed as part of a thesis and solves for both the steady-state and transition dynamics of agents across multiple locations.

## Repository Structure

- `Main.m` — Entry point that sets parameters, solves the no-help equilibrium, evaluates the benchmark and policy counterfactuals, and computes welfare/decomposition statistics.
- `core/` — Functions that implement the dynamic programming routines (e.g. `PolicyDynamics`, `TransitionDynamics`, `noHelpEqm`).
- `utils/` — Helper routines for parameter initialization, grid creation, interpolation, and agent simulation.

## Requirements

- MATLAB R2023b or later.
- The code uses only base MATLAB functions; no additional toolboxes should be required.

## Getting Started

1. Clone this repository and open `Main.m` in MATLAB.
2. Adjust the welfare horizon `analysisHorizon` (\tilde{T}) in Section 4 if needed. By default it is set to 40 periods, which is shorter than the simulation horizon `settings.T`.
3. Run `Main.m` to solve each scenario (benchmark plus transport and shelter programs with different start dates) and to produce the welfare comparisons.

Running `Main.m` prints progress messages, saves location-share plots in the `results/` folder, and writes welfare summaries to `results/welfare/welfare_summary.csv`. A `.mat` file with the detailed scenario outputs is saved in the same directory. The default parameters and iteration settings—including the random seed used for matching—can be adjusted by editing the files in the `utils/` directory.

## Citation

If you use this code in your own work, please cite the thesis or contact the author. The functions include detailed comments and were written by Agustin Deambrosi.
