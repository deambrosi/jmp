function matrices = constructMatrix(dims, params, grids, indexes)
% CONSTRUCTMATRIX Precomputes matrices for utility and after-migration wealth.
%
%   INPUTS:
%       dims        - Struct with dimension parameters (N, K, etc.)
%       params      - Struct with model parameters (bbeta, A, B, ssigma, ttau)
%       grids       - Struct with state grids (agrid, ahgrid, eta, psi)
%       indexes     - Struct with index matrices for reshaping and grid use
%
%   OUTPUT:
%       matrices    - Struct containing:
%                       .Ue       - Utility from consumption while employed
%                       .a_prime  - Coarse grid asset levels net of migration costs
%                       .A_prime  - Asset post-migration grid expanded across states
%
%   AUTHOR: Agustin Deambrosi
%   LAST REVISED: April 2025
% =========================================================================

    %% 1. Compute Consumption and Utility (Ue)
    cons			= (1 / params.bbeta) .* grids.agrid(indexes.I_ap) + ...
                      params.A(indexes.I_Np) .* grids.eta(indexes.I_ep) - ...
                      grids.ahgrid(indexes.I_app);

    amenity_weight	= params.B(indexes.I_Np) .* grids.psi(indexes.I_psip);
    %Uu				= log(cons());
    
    Ue              =   zeros(size(cons));
    Ue(cons > 0)	=   amenity_weight(cons > 0) .* log(cons(cons > 0));
    Ue(cons <= 0)	=   -realmax;   % Penalize infeasible consumption
                    
    %% 2. Compute After-Migration Wealth Matrix (a' and A')
    mig_costs		= permute(params.ttau, [4, 1, 2, 3]);      % Adds singleton 4th dim: wealth
    a_prime			= grids.agrid - mig_costs;                  % Adjust for migration cost

    A_Prime			= permute(a_prime, [5, 1, 2, 3, 4]);       % Expand to 5th Dim: 'state'
    A_Prime			= repmat(A_Prime, [dims.S, 1, 1, 1, 1]);   % Copy across productivity states

    %% 3. Output Struct
    matrices.Ue			= Ue;
    matrices.a_prime	= a_prime;
    matrices.A_prime	= A_Prime;

end
