function [vf_path, pol_path] = PolicyDynamicsVarCost(M1, vf_terminal, dims, params, grids, indexes, matrices_path, settings)
% POLICYDYNAMICSVARCOST Backward induction with time-varying migration costs.
%
%   Inputs are the same as PolicyDynamics, with the addition of
%   'matrices_path', a cell array of length T containing the precomputed
%   matrices for each period.
%
%   AUTHOR: ChatGPT
%   DATE: 2025
% =========================================================================

    T                  = settings.T;
    vf_path            = cell(T, 1);
    pol_path.a         = cell(T, 1);
    pol_path.an        = cell(T, 1);
    pol_path.mu        = cell(T, 1);
    pol_path.mun       = cell(T, 1);

    G_path0            = computeG(M1, params.ggamma);
    vf_path{T}         = vf_terminal;

    for t = T-1:-1:1
        matrices_t      = matrices_path{t};
        G_path_t        = G_path0(:, t);

        [vf_t, pol_t]   = updateValueAndPolicy(vf_path{t+1}, dims, params, grids, indexes, matrices_t, G_path_t);

        vf_path{t}       = vf_t;
        pol_path.a{t}    = pol_t.a;
        pol_path.an{t}   = pol_t.an;
        pol_path.mu{t}   = pol_t.mu;
        pol_path.mun{t}  = pol_t.mun;
    end
end
