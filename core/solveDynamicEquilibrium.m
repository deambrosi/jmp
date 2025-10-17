function [pol_eqm, M_eqm, it_count, vf_path] = solveDynamicEquilibrium(M0, vf_terminal, m0, dims, params, grids, indexes, matrices, settings, aidParams)
% SOLVEDYNAMICEQUILIBRIUM Solves for the dynamic migration equilibrium.
%
%   Iterates over the dynamic transition path of the economy using:
%       - Backward induction for value and policy functions
%       - Forward simulation to update agent distributions
%       - Weighted convergence criterion favoring early periods
%
%   INPUTS:
%       M0           - [N x T] initial guess for network agent distribution path
%       vf_terminal  - Struct with steady-state value functions (.V, .Vn, .R, .Rn)
%       m0           - Initial agent distribution (struct array)
%       dims         - Model dimensions (struct)
%       params       - Model parameters (struct)
%       grids        - Grids (struct)
%       indexes      - Indexing structures (struct)
%       matrices     - Precomputed matrix structures (Ue, a_prime, A_prime)
%       settings     - Simulation and iteration settings (struct)
%       aidParams    - (optional) struct to apply temporary aid when computing
%                      the forward simulation. If empty or omitted, no aid is
%                      used.
%
%   OUTPUTS:
%       pol_eqm   - Struct with converged dynamic policy functions (.a, .an, .mu, .mun)
%       M_eqm     - [N x T] converged path of network agent distribution
%       it_count  - Number of iterations until convergence
%       vf_path   - Cell array of value functions over time (from final iteration)
%
%   AUTHOR: Agustin Deambrosi
%   LAST REVISED: April 2025
% =========================================================================

    %% 1. Initialization
    T               = settings.T;
    diffM           = 1;
    it_count        = 0;
    M_eqm           = M0;

    % Exponential weights: early periods get higher weight
    beta_weight     = 0.7;                             % Adjust as needed
    time_weights    = beta_weight .^ (0:T-1);          % [1 x T]
    time_weights    = time_weights / sum(time_weights);% Normalize to sum to 1

    %% 2. Iteration Loop
    fprintf('\nSolving Dynamic Equilibrium...\n');
    while (diffM > settings.tolM) && (it_count < settings.MaxItJ)
        it_count = it_count + 1;

        % Step 1: Policy functions via backward induction
        [vf_path, pol_new] = PolicyDynamics(M_eqm, vf_terminal, dims, params, grids, indexes, matrices, settings);

        % Step 2: Forward simulation of agent dynamics
        if nargin < 10 || isempty(aidParams)
            [~, M_new, ~] = simulateAgents(m0, pol_new, repmat(params.G0, 1, T), ...
                                           dims, params, grids, settings);
        else
            [~, M_new, ~, ~] = simulateAgentsWithAid(m0, pol_new, ...
                repmat(params.G0, 1, T), dims, params, grids, settings, aidParams);
        end

        % Step 3: Compute weighted difference across periods
        diff_per_t = sum(abs(M_eqm - M_new) ./ (1 + abs(M_new)), 1);   % [1 x T]
        diffM      = sum(time_weights .* diff_per_t);                 % Scalar

        % Step 4: Update
        M_eqm      = M_new;

        % Step 5: Print status
        fprintf('  Iteration %d: Weighted diffM = %.6f\n', it_count, diffM);
    end

    % Final output
    pol_eqm = pol_new;
    % vf_path already contains value functions from the last iteration
end
