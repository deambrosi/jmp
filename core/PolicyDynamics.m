function [vf_path, pol_path] = PolicyDynamics(M1, vf_terminal, dims, params, grids, indexes, matrices, settings, helpOverride)
% POLICYDYNAMICS Computes policy functions over the transition path.
%
%   [vf_path, pol_path] = PolicyDynamics(M1, vf_terminal, ..., helpOverride)
%
%   INPUTS:
%       M1         - [N x T] matrix with guessed path of network agent distribution
%       vf_terminal- Struct with steady-state value functions at T (.V, .Vn, etc.)
%       dims       - Model dimensions (struct)
%       params     - Model parameters (includes .ggamma, .ttau, .cchi, etc.)
%       grids      - Grid structure (.agrid, .ahgrid, etc.)
%       indexes    - Index structure for looping and vectorization
%       matrices   - Struct with Ue, a_prime, A_prime
%       settings   - Simulation/iteration controls
%
%   OPTIONAL INPUT:
%       helpOverride - [H x T] matrix supplying exogenous help distributions
%                      to be used at each date. When provided, the path is
%                      used verbatim instead of computing it from M1. This is
%                      useful for shutting down network effects while keeping
%                      the rest of the equilibrium logic intact.
%
%   OUTPUTS:
%       vf_path    - Cell array of value function structs over time
%       pol_path   - Struct with dynamic policies:
%                     .a, .an, .mu, .mun for each t
%
%   AUTHOR: Agustin Deambrosi
%   LAST REVISED: April 2025
% =========================================================================

    %% 1. Setup
    T                  = settings.T;
    vf_path            = cell(T, 1);
    pol_path.a         = cell(T, 1);
    pol_path.an        = cell(T, 1);
    pol_path.mu        = cell(T, 1);
    pol_path.mun       = cell(T, 1);

    if nargin < 9 || isempty(helpOverride)
        G_path0 = computeG(M1, params.ggamma);   % help vector distributions from masses
    else
        G_path0 = helpOverride;
        if size(G_path0, 2) ~= T
            error('helpOverride must have T = %d columns.', T);
        end
    end

    vf_path{T}         = vf_terminal;  % Terminal condition

    %% 2. Backward Induction: t = T-1 to 1
    for t = T-1:-1:1
        % A) Compute help probabilities G_t based on M1(:, t)
        G_path_t      = G_path0(:,t);  % [H x 1]

        % B) Use continuation value at t+1 to update policy and value function at t
        [vf_t, pol_t]  = updateValueAndPolicy( ...
            vf_path{t+1}, dims, params, grids, indexes, matrices, G_path_t);

        % C) Store results
        vf_path{t}       = vf_t;
        pol_path.a{t}    = pol_t.a;
        pol_path.an{t}   = pol_t.an;
        pol_path.mu{t}   = pol_t.mu;
        pol_path.mun{t}  = pol_t.mun;
    end
end
