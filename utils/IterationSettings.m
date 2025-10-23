function settings = IterationSettings()
% ITERATIONSETTINGS Initializes iteration parameters for solving the model.
%
%   OUTPUT:
%       settings - Struct containing:
%           .it       - Iteration counter (initially 0)
%           .diffV    - Initial V-function difference
%           .tolV     - Convergence tolerance for V
%           .tolM     - Tolerance for migration convergence
%           .MaxItV   - Max iterations for value function iteration
%           .MaxItJ   - Max iterations for the outer equilibrium loop
%           .MaxIter  - Safety cap on miscellaneous algorithms
%           .Nagents  - Number of simulated agents
%           .T        - Total simulation horizon
%           .burn     - Burn-in periods removed from moment calculations
%           .rngSeed  - Seed used to synchronize the Monte Carlo draws
%
%   AUTHOR: Agustin Deambrosi (updated by ChatGPT)
%   LAST REVISED: May 2025
% =========================================================================

    %% Iteration counters and convergence values
    settings.it         = 0;
    settings.diffV      = 1;

    %% Tolerances
    settings.tolV       = 0.5;
    settings.tolM       = 1e-2;

    %% Iteration limits
    settings.MaxItV     = 40;
    settings.MaxItJ     = 2;
    settings.MaxIter    = 100;

    %% Simulation settings
    settings.Nagents    = 5000;
    settings.T          = 100;
    settings.burn       = 50;

    %% Random-number generator seed for reproducibility across scenarios
    settings.rngSeed    = 12345;

end

