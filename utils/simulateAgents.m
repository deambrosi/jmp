function [M_history, MIN_history, agentData] = simulateAgents(m0, pol, G_dist, dims, params, grids, settings)
% SIMULATEAGENTS Simulates agent evolution over T periods using policy functions.
%
%   [M_history, MIN_history, agentData] = simulateAgents(...)
%
%   Simulates agent dynamics over T periods:
%       - Applies policy functions for saving and migration.
%       - Evolves productivity and network status.
%       - Tracks all relevant state variables for each agent.
%
%   INPUTS:
%       m0        - [Nagents x 1] initial agent structs with fields:
%                     .state, .wealth, .location, .network
%       pol       - Policy struct with fields:
%                     .a, .an, .mu, .mun (each T-1 cell array)
%       T         - Simulation horizon (number of periods)
%       G_dist    - [H x T] matrix with time-varying help vector PMF
%       dims      - Dimension struct (N, Na, S, H)
%       params    - Parameter struct (.ttau, .P, .cchi)
%       grids     - Struct with asset grids (.agrid, .ahgrid)
%       settings  - Struct with settings (.Nagents)
%
%   OUTPUTS:
%       M_history   - [N x T] total agent count per location per period
%       MIN_history - [N x T] networked agent count per location per period
%       agentData   - Struct with full trajectories:
%                        .location, .wealth, .state, .network
%
%   AUTHOR: Agustin Deambrosi
%   LAST REVISED: April 2025
% =========================================================================

    %% 1. Setup
    T   =   settings.T;
    
    isTimeInvariant = ~iscell(pol.a);
    if isTimeInvariant
        pol.a   = repmat({pol.a}, T-1, 1);
        pol.an  = repmat({pol.an}, T-1, 1);
        pol.mu  = repmat({pol.mu}, T-1, 1);
        pol.mun = repmat({pol.mun}, T-1, 1);
    end

    numAgents    = settings.Nagents;
    numLocations = dims.N;

    %% 2. Preallocate full agent trajectories
    locationTraj = zeros(numAgents, T);
    wealthTraj   = zeros(numAgents, T);
    stateTraj    = zeros(numAgents, T);
    networkTraj  = zeros(numAgents, T);

    %% 3. Simulation Loop
    parfor agentIdx = 1:numAgents
        % Initialize agent
        agent = m0(agentIdx);
        locHist = zeros(1, T);
        weaHist = zeros(1, T);
        staHist = zeros(1, T);
        netHist = zeros(1, T);

        % Initial conditions
        loc  = agent.location;
        wea  = agent.wealth;
        sta  = agent.state;
        net  = agent.network;

        locHist(1) = loc;
        weaHist(1) = wea;
        staHist(1) = sta;
        netHist(1) = net;

        for t = 1:(T-1)
            % A) Asset decision
            if net == 1
                a_fine = pol.an{t}(sta, wea, loc);
            else
                a_fine = pol.a{t}(sta, wea, loc);
            end
            a_fine = min(max(a_fine, 1), length(grids.ahgrid));
            [~, nextWea] = min(abs(grids.agrid - grids.ahgrid(a_fine)));

            % B) Migration decision
            if net == 1
                G_t     = G_dist(:, t);
                h_idx   = find(cumsum(G_t) >= rand(), 1);
                migProb = squeeze(pol.mun{t}(sta, wea, loc, :, h_idx));
            else
                migProb = squeeze(pol.mu{t}(sta, wea, loc, :));
                h_idx   = 1;  % default help vector
            end
            migProb = migProb / sum(migProb);
            nextLoc = find(cumsum(migProb) >= rand(), 1);
            
            % fallback if find returns empty (due to numerical issues)
            if isempty(nextLoc)
                [~, nextLoc] = max(migProb);
            end

            % C) Wealth and state transitions
            if nextLoc ~= loc
                migCost     = params.ttau(loc, nextLoc, h_idx);
                newA        = grids.agrid(nextWea) - migCost;
                [~, nextWea] = min(abs(grids.agrid - newA));
                sta         = 1;  % reset (η, ψ) after migration
            else
                transCdf    = cumsum(params.P(sta, :));
                sta         = find(transCdf >= rand(), 1);
            end

            % D) Network status update
            if net == 1 & nextLoc ~= 1
                if rand() < params.cchi
                    net = 0;
                end
            end

            % Save next-period state
            loc = nextLoc;
            wea = nextWea;

            locHist(t+1) = loc;
            weaHist(t+1) = wea;
            staHist(t+1) = sta;
            netHist(t+1) = net;
        end

        % Save trajectories
        locationTraj(agentIdx, :) = locHist;
        wealthTraj(agentIdx, :)   = weaHist;
        stateTraj(agentIdx, :)    = staHist;
        networkTraj(agentIdx, :)  = netHist;
    end

    %% 4. Aggregate Location Histories
    M_history = zeros(numLocations, T);
    MIN_history = zeros(numLocations, T);

    for t = 1:T
        locs = locationTraj(:, t);
        nets = networkTraj(:, t);

        M_history(:, t) = accumarray(locs, 1, [numLocations, 1]) / numAgents;
        MIN_history(:, t) = accumarray(locs(nets == 1), 1, [numLocations, 1]) / numAgents;
    end

    %% 5. Output full agent trajectories
    agentData.location = locationTraj;
    agentData.wealth   = wealthTraj;
    agentData.state    = stateTraj;
    agentData.network  = networkTraj;
end
