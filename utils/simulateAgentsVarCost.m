function [M_history, MIN_history, agentData] = simulateAgentsVarCost(m0, pol, G_dist, dims, params, grids, settings, tau_path)
% SIMULATEAGENTSVARCOST Simulates agent dynamics with time-varying migration costs.
%
%   This is an extension of simulateAgents where migration costs can vary
%   by period. tau_path is a [N x N x H x T] tensor specifying the cost for
%   each period t.
%
%   All other inputs are identical to simulateAgents.
%
%   AUTHOR: ChatGPT
%   DATE: 2025
% =========================================================================

    T   = settings.T;

    isTimeInvariant = ~iscell(pol.a);
    if isTimeInvariant
        pol.a   = repmat({pol.a}, T-1, 1);
        pol.an  = repmat({pol.an}, T-1, 1);
        pol.mu  = repmat({pol.mu}, T-1, 1);
        pol.mun = repmat({pol.mun}, T-1, 1);
    end

    numAgents    = settings.Nagents;
    numLocations = dims.N;

    locationTraj = zeros(numAgents, T);
    wealthTraj   = zeros(numAgents, T);
    stateTraj    = zeros(numAgents, T);
    networkTraj  = zeros(numAgents, T);

    parfor agentIdx = 1:numAgents
        agent = m0(agentIdx);
        locHist = zeros(1, T);
        weaHist = zeros(1, T);
        staHist = zeros(1, T);
        netHist = zeros(1, T);

        loc  = agent.location;
        wea  = agent.wealth;
        sta  = agent.state;
        net  = agent.network;

        locHist(1) = loc;
        weaHist(1) = wea;
        staHist(1) = sta;
        netHist(1) = net;

        for t = 1:(T-1)
            if net == 1
                a_fine = pol.an{t}(sta, wea, loc);
            else
                a_fine = pol.a{t}(sta, wea, loc);
            end
            a_fine = min(max(a_fine, 1), length(grids.ahgrid));
            [~, nextWea] = min(abs(grids.agrid - grids.ahgrid(a_fine)));

            if net == 1
                G_t     = G_dist(:, t);
                h_idx   = find(cumsum(G_t) >= rand(), 1);
                migProb = squeeze(pol.mun{t}(sta, wea, loc, :, h_idx));
            else
                migProb = squeeze(pol.mu{t}(sta, wea, loc, :));
                h_idx   = 1;
            end
            migProb = migProb / sum(migProb);
            nextLoc = find(cumsum(migProb) >= rand(), 1);
            if isempty(nextLoc)
                [~, nextLoc] = max(migProb);
            end

            if nextLoc ~= loc
                migCost     = tau_path(loc, nextLoc, h_idx, t);
                newA        = grids.agrid(nextWea) - migCost;
                [~, nextWea] = min(abs(grids.agrid - newA));
                sta         = 1;
            else
                transCdf    = cumsum(params.P(sta, :));
                sta         = find(transCdf >= rand(), 1);
            end

            if net == 1 && nextLoc ~= 1
                if rand() < params.cchi
                    net = 0;
                end
            end

            loc = nextLoc;
            wea = nextWea;

            locHist(t+1) = loc;
            weaHist(t+1) = wea;
            staHist(t+1) = sta;
            netHist(t+1) = net;
        end

        locationTraj(agentIdx, :) = locHist;
        wealthTraj(agentIdx, :)   = weaHist;
        stateTraj(agentIdx, :)    = staHist;
        networkTraj(agentIdx, :)  = netHist;
    end

    M_history = zeros(numLocations, T);
    MIN_history = zeros(numLocations, T);
    for t = 1:T
        locs = locationTraj(:, t);
        nets = networkTraj(:, t);
        M_history(:, t) = accumarray(locs, 1, [numLocations, 1]) / numAgents;
        MIN_history(:, t) = accumarray(locs(nets == 1), 1, [numLocations, 1]) / numAgents;
    end

    agentData.location = locationTraj;
    agentData.wealth   = wealthTraj;
    agentData.state    = stateTraj;
    agentData.network  = networkTraj;
end
