function m0 = createInitialDistribution(dims, settings)
% CREATEINITIALDISTRIBUTION Generates initial agent distribution over states.
%
%   OUTPUT:
%       m0 - [Nagents x 1] struct array with fields:
%               .state     - Joint state index s ∈ {1,...,dims.S} for (η, ψ)
%               .wealth    - Wealth index on coarse grid (1 to dims.Na)
%               .location  - Initial location (all agents start in 1)
%               .network   - Network affiliation: 1 = in network, 0 = out
%
%   INPUTS:
%       dims     - Struct with model dimensions (S, Na, N)
%       settings - Struct with simulation control (Nagents: number of agents)
%
%   AUTHOR: Agustin Deambrosi
%   LAST REVISED: April 2025
% =========================================================================

    numAgents	= settings.Nagents;

    % Preallocate agent structure
    m0			= repmat(struct('state', 0, 'wealth', 0, 'location', 0, 'network', 0), numAgents, 1);

    % Initialize agents randomly over state and wealth; fixed location and network status
    for i = 1:numAgents
        m0(i).state		= randi(dims.S);      % Joint state (η, ψ)
        m0(i).wealth	= randi(dims.Na);     % Asset index on coarse grid
        m0(i).location	= 1;                  % Start in location 1
        m0(i).network	= 1;                  % All start in the network
    end

end
