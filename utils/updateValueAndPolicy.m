function [vf, pol] = updateValueAndPolicy(val, dims, params, grids, indexes, matrices, G)
% UPDATEVALUEANDPOLICY Updates value functions and optimal asset/migration policies.
%
%   This function computes updated value functions for agents who:
%       - Are **not** affiliated with the migrant network (n = 0)
%       - Are still in the network and may receive help offers (n = 1)
%
%   Steps:
%       1. Compute continuation values if agents stay in place.
%       2. Compute continuation values if agents migrate (with interpolation).
%       3. Calculate softmax-based migration probabilities over destinations.
%       4. Take expectations over future states and help offers.
%       5. Compute the expected value functions (R, Rn).
%       6. Find optimal asset choices via maximization.
%       7. Construct migration policy functions (mu, mun).
%
%
%   INPUTS:
%       val      - Struct with current value functions: .V  (n = 0), .Vn (n = 1)
%       dims     - Model dimensions (S, Na, N, H, etc.)
%       params   - Model parameters (bbeta, cchi, CONS, nnu, etc.)
%       grids    - Grid variables (agrid, ahgrid, etc.)
%       indexes  - Index matrices (I_ep, I_ap, II, etc.)
%       matrices - Struct with utility matrix and after-migration wealth: Ue, a_prime, A_prime
%       G        - Help offer PMF over h ∈ {0,1}^N (only relevant for n = 1)
%
%   OUTPUTS:
%       vf       - Struct with updated value functions:
%                    .V   : Value function for non-network agents
%                    .Vn  : Value function for network agents
%                    .R   : Continuation value for non-network
%                    .Rn  : Continuation value for network
%
%       pol      - Struct with optimal policy rules:
%                    .a   : Asset choice index (n = 0)
%                    .an  : Asset choice index (n = 1)
%                    .mu  : Migration prob. to each destination (n = 0)
%                    .mun : Migration prob. to each destination (n = 1)
%
%   AUTHOR: Agustin Deambrosi
%   LAST REVISED: April 2025
% =========================================================================

    %% 1. Continuation Values — if agent stays in place
    cont_no_mig			= pagemtimes(params.P, val.V);  % n = 0
    cont_no_mig_net		= pagemtimes(params.P, ...
                          (1 - params.cchi) * val.Vn + params.cchi * val.V); % n = 1

    %% 2. Continuation Values — if agent migrates
    cont_mig			= interp_migration(grids.agrid, matrices.a_prime, dims, val.V);
    cont_mig_net		= interp_migration(grids.agrid, matrices.a_prime, dims, ...
                             (1 - params.cchi) * val.Vn + params.cchi * val.V);

    % Assign value of staying (i → i) across all help configurations h
    cont_mig(indexes.II)		= repmat(cont_no_mig, 1, 1, 1, dims.H);
    cont_mig_net(indexes.II)	= repmat(cont_no_mig_net, 1, 1, 1, dims.H);

    % Make infeasible migration moves very unattractive
    cont_mig(matrices.A_prime < 0)		= -Inf;
    cont_mig_net(matrices.A_prime < 0)	= -Inf;

    %% 3. Migration Probabilities (softmax over destinations)
    exp_vals		= (exp(cont_mig / params.CONS)).^(1 / params.nnu);
    exp_vals_net	= (exp(cont_mig_net / params.CONS)).^(1 / params.nnu);

    mmuu			= exp_vals ./ sum(exp_vals, 4);        % n = 0
    mmuu_net		= exp_vals_net ./ sum(exp_vals_net, 4);% n = 1

    %% 4. Expected Continuation Value
    % For agents not in network, only one possible help vector (all zeros)
    exp_value					= mmuu(:,:,:,:,1) .* cont_mig(:,:,:,:,1);
    exp_value(isnan(exp_value)) = 0;           

    % For network agents: weighted average over help vectors h
    exp_value_net				= mmuu_net .* cont_mig_net;
    exp_value_net(isnan(exp_value_net)) = 0;                % Handle numerical artifacts
    exp_value_net				= permute(exp_value_net, [5, 1, 2, 3, 4]);    % h → first dim
    exp_value_net				= sum(exp_value_net .* reshape(G, [dims.H 1 1 1 1]), 1);
    exp_value_net				= permute(exp_value_net, [2, 3, 4, 5, 1]);    % h → last dim

    %% 5. Compute Expected Value: R = β E[V']
    vf.R						= params.bbeta * sum(exp_value, 4);         % n = 0
    vf.Rn						= params.bbeta * sum(exp_value_net, 4);     % n = 1

    %% 6. Asset Choice: maximize utility + continuation value
    interp_R					= interpolateToFinerGrid(grids.agrid, grids.ahgrid, vf.R);
    total_val					= matrices.Ue + interp_R;
    [vf.V, pol.a]				= max(total_val, [], 4);   % n = 0

    interp_Rn					= interpolateToFinerGrid(grids.agrid, grids.ahgrid, vf.Rn);
    total_valn					= matrices.Ue + interp_Rn;
    [vf.Vn, pol.an]				= max(total_valn, [], 4);  % n = 1

    %% 7. Migration Flows
    pol.mu						= mmuu(:,:,:,:,1);    % Migration policy (n = 0): only one help vector
    pol.mun						= mmuu_net;

end
