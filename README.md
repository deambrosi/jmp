# Dynamic Migration Network Model

This repository contains MATLAB code for simulating a dynamic migration model with network-based help. The code was developed as part of a thesis and solves for both the steady-state and transition dynamics of agents across multiple locations.

## Repository Structure

- `Main.m` — Entry point that sets parameters, solves the no-help equilibrium, computes the dynamic equilibrium, and runs final simulations.
- `core/` — Functions that implement the dynamic programming routines (e.g. `PolicyDynamics`, `TransitionDynamics`, `noHelpEqm`).
- `utils/` — Helper routines for parameter initialization, grid creation, interpolation, and agent simulation.

## Requirements

- MATLAB R2023b or later.
- The code uses only base MATLAB functions; no additional toolboxes should be required.

## Getting Started

1. Clone this repository and open `Main.m` in MATLAB.
2. Run `Main.m` to initialize parameters and compute the steady-state and transition dynamics.

Running `Main.m` prints progress messages and generates simulation results. The default parameters and iteration settings can be adjusted by editing the files in the `utils/` directory.

## Citation

If you use this code in your own work, please cite the thesis or contact the author. The functions include detailed comments and were written by Agustin Deambrosi.
