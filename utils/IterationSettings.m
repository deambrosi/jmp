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
%           .MaxIter  - Max iterations for outer algorithm
%           .Nagents  - Number of agents in simulation
%           .T        - Total time periods simulated
%           .burn     - Burn-in periods removed before analysis
%
%   AUTHOR: Agustin Deambrosi
%   LAST REVISED: April 2025
% =========================================================================

    %% Iteration counters and convergence values
    settings.it         = 0;
    settings.diffV      = 1;

    %% Tolerances
    settings.tolV       = 0.5;
    settings.tolM       = 1e-2;

    %% Iteration limits
    settings.MaxItV     = 40;
    settings.MaxItJ     = 3;
    settings.MaxIter    = 100;

    %% Simulation settings
    settings.Nagents    = 5000;
    settings.T          = 100;
    settings.burn       = 50;

    %% Display controls
    settings.showDynamicEqmProgress = false;

end
