function [pol_eqm, M_eqm, it_count, vf_path] = solveDynamicEquilibriumVarCost(M0, vf_terminal, m0, dims, params, grids, indexes, matrices_path, tau_path, settings)
% SOLVEDYNAMICEQUILIBRIUMVARCOST Dynamic equilibrium with time-varying migration costs.
%
%   This mirrors solveDynamicEquilibrium but allows migration costs to vary
%   over time. matrices_path and tau_path should be constructed using
%   createMigrationCostPath.
%
%   AUTHOR: ChatGPT
%   DATE: 2025
% =========================================================================

    T               = settings.T;
    diffM           = 1;
    it_count        = 0;
    M_eqm           = M0;

    beta_weight     = 0.7;
    time_weights    = beta_weight .^ (0:T-1);
    time_weights    = time_weights / sum(time_weights);

    fprintf('\nSolving Dynamic Equilibrium with migration cost shock...\n');
    while (diffM > settings.tolM) && (it_count < settings.MaxItJ)
        it_count = it_count + 1;

        [vf_path, pol_new] = PolicyDynamicsVarCost(M_eqm, vf_terminal, dims, params, grids, indexes, matrices_path, settings);

        [~, M_new, ~] = simulateAgentsVarCost(m0, pol_new, repmat(params.G0, 1, T), dims, params, grids, settings, tau_path);

        diff_per_t = sum(abs(M_eqm - M_new) ./ (1 + abs(M_new)), 1);
        diffM      = sum(time_weights .* diff_per_t);
        M_eqm      = M_new;

        fprintf('  Iteration %d: Weighted diffM = %.6f\n', it_count, diffM);
    end

    pol_eqm = pol_new;
end
