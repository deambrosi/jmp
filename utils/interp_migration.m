function f = interp_migration(grid, xgrid, dims, V)
% INTERP_MIGRATION Interpolates value function across post-migration asset levels.
%
%   This function takes a value function V defined over (eta, a, i), and 
%   interpolates it over asset values a - tau for each possible destination
%   and help vector h. This is used to compute continuation values after
%   migration decisions are made.
%
%   INPUTS:
%       grid   - [Na x 1] Coarse wealth grid
%       xgrid  - [Na x N x N x H] Grid of post-migration asset levels (a - tau)
%       dims   - Struct with model dimensions (Na, N, H, S, etc.)
%       V      - [K x Na x N] Value function before migration
%
%   OUTPUT:
%       f      - [S x Na x N x N x H] Interpolated continuation values
%
%   AUTHOR: Agustin Deambrosi
%   LAST REVISED: April 2025
% =========================================================================

    %% 1. Rearrange and reshape V for interpolation
    % Move asset to first dimension: [Na x N x K]
    V_reordered		= permute(V, [2, 3, 1]); 

    % We'll interpolate based on productivity level 1 (assumed representative)
    % This gives us: [Na x N]
    V_base			= V_reordered(:, :, 1); 

    % Expand to [Na x current_loc x destination_loc x help_vector]
    V_expanded		= repmat(V_base, 1, 1, dims.N, dims.H);
    V_expanded		= permute(V_expanded, [1, 3, 2, 4]);

    % Flatten for batch interpolation
    V_flat			= reshape(V_expanded, dims.Na, []);
    xgrid_flat		= reshape(xgrid, dims.Na, []);

    %% 2. Interpolate in parallel over all destination/help configurations
    f_interp		= zeros(size(V_flat));
    parfor j = 1:size(V_flat, 2)
        f_interp(:, j) = interp1(grid, V_flat(:, j), xgrid_flat(:, j), 'linear', 'extrap');
    end

    %% 3. Reshape to [Na x N x N x H] and expand across S (state space)
    f_shaped		= reshape(f_interp, dims.Na, dims.N, dims.N, dims.H);
    f_shaped		= permute(f_shaped, [5, 1, 2, 3, 4]);  % add singleton dimension for state

    % Replicate across S to match full dimensionality [S x Na x N x N x H]
    f				= repmat(f_shaped, [dims.S, 1, 1, 1, 1]);

end
