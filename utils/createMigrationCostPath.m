function [tau_path, matrices_path] = createMigrationCostPath(params, dims, grids, indexes, settings, factor, t_start, duration)
% CREATEMIGRATIONCOSTPATH Construct migration cost and matrix paths with a temporary shock.
%
%   [tau_path, matrices_path] = createMigrationCostPath(params, dims, grids, indexes, settings, factor, t_start, duration)
%
%   Constructs a time-varying migration cost tensor and a cell array of
%   precomputed matrices for dynamic programming. The migration costs are
%   reduced by 'factor' starting at period t_start for 'duration' periods.
%
%   INPUTS:
%       params      - Parameter struct with baseline params.ttau (after help)
%       dims        - Dimension struct
%       grids       - Grids for constructMatrix
%       indexes     - Index structures for constructMatrix
%       settings    - Settings struct with .T (number of periods)
%       factor      - Multiplicative factor applied during the shock
%       t_start     - First period of the shock (1-indexed)
%       duration    - Number of periods the shock lasts
%
%   OUTPUTS:
%       tau_path       - [N x N x H x T] tensor of migration costs per period
%       matrices_path  - 1xT cell array with precomputed matrices per period
%
%   AUTHOR: ChatGPT
%   DATE: 2025
% =========================================================================

    T           = settings.T;
    base_tau    = params.ttau;

    % Precompute baseline matrix using existing parameters
    matrices_base   = constructMatrix(dims, params, grids, indexes);

    % Construct low cost parameters and matrices
    params_low      = params;
    params_low.ttau = factor .* base_tau;
    matrices_low    = constructMatrix(dims, params_low, grids, indexes);

    % Initialize outputs
    tau_path       = repmat(base_tau, 1, 1, 1, T);
    matrices_path  = repmat({matrices_base}, T, 1);

    t_end = min(t_start + duration - 1, T);
    shock_len = t_end - t_start + 1;
    tau_path(:,:,:,t_start:t_end) = repmat(factor .* base_tau, 1, 1, 1, shock_len);
    for tt = t_start:t_end
        matrices_path{tt} = matrices_low;
    end
end
