function decomposition = decomposeWelfareDifference(counter, baseline, beta, matchIdx, timingCounter, timingBaseline, varargin)
% DECOMPOSEWELFAREDIFFERENCE Implements the three-way welfare decomposition.
%
%   decomposition = decomposeWelfareDifference(counter, baseline, beta,
%       matchIdx, timingCounter, timingBaseline, welfareHorizon) returns agent-
%   level and aggregate contributions for the welfare difference W^c - W^b
%   following the formulas provided in the project brief.
%
%   INPUTS:
%       counter, baseline - Structs from computeScenarioWelfare containing
%                           .valuePanel and .betaPowers
%       beta              - Discount factor β
%       matchIdx          - Mapping f(i) pairing counterfactual agents with
%                           benchmark agents
%       timingCounter     - Struct from computeMigrationTiming for scenario c
%       timingBaseline    - Struct from computeMigrationTiming for scenario b
%       welfareHorizon    - (Optional) Horizon \tilde{T}. If omitted the
%                           function falls back to the welfare horizon stored
%                           in the welfare structs.
%
%   OUTPUT:
%       decomposition - Struct with fields:
%                           .perAgent.leaveEarlier
%                           .perAgent.pathQuality
%                           .perAgent.destination
%                           .perAgent.totalDifference
%                           .aggregate (same components summed over agents)
%
%   AUTHOR: ChatGPT
%   DATE:   May 2025
% =========================================================================

    numAgents = size(counter.valuePanel, 1);

    if ~isempty(varargin)
        Ttilde = varargin{1};
    elseif isfield(counter, 'horizon')
        Ttilde = counter.horizon;
    else
        Ttilde = size(counter.valuePanel, 2);
    end

    if size(counter.valuePanel, 2) < Ttilde || size(baseline.valuePanel, 2) < Ttilde
        error('Requested welfare horizon exceeds available value data.');
    end

    betaVec = beta .^ (1:Ttilde);

    perAgent.leaveEarlier   = zeros(numAgents, 1);
    perAgent.pathQuality    = zeros(numAgents, 1);
    perAgent.destination    = zeros(numAgents, 1);
    perAgent.totalDifference= zeros(numAgents, 1);

    for iAgent = 1:numAgents
        jAgent = matchIdx(iAgent);

        vc = counter.valuePanel(iAgent, 1:Ttilde);
        vb = baseline.valuePanel(jAgent, 1:Ttilde);

        diffSeries = vc - vb;
        perAgent.totalDifference(iAgent) = sum(betaVec .* diffSeries);

        if ~(timingCounter.everLeft(iAgent) && timingBaseline.everLeft(jAgent))
            perAgent.destination(iAgent) = perAgent.totalDifference(iAgent);
            continue;
        end

        deltaTstar = timingBaseline.tStar(jAgent) - timingCounter.tStar(iAgent);
        vbShifted  = shiftSeriesWithPadding(vb, deltaTstar);
        leaveTerm  = vbShifted - vb;
        perAgent.leaveEarlier(iAgent) = sum(betaVec .* leaveTerm);

        gapSeries = vc - vbShifted;

        tMaxDouble = max(timingBaseline.tDoubleStar(jAgent), timingCounter.tDoubleStar(iAgent));
        cutoff     = min(Ttilde, timingCounter.tStar(iAgent) + tMaxDouble);

        earlyIdx = 1:cutoff;
        lateIdx  = (cutoff + 1):Ttilde;

        if ~isempty(earlyIdx)
            perAgent.pathQuality(iAgent) = sum(betaVec(earlyIdx) .* gapSeries(earlyIdx));
        end

        if ~isempty(lateIdx)
            perAgent.destination(iAgent) = sum(betaVec(lateIdx) .* gapSeries(lateIdx));
        end

        % If cutoff == Ttilde the destination component should receive zero,
        % which is already achieved because lateIdx is empty.
    end

    aggregate.leaveEarlier    = sum(perAgent.leaveEarlier);
    aggregate.pathQuality     = sum(perAgent.pathQuality);
    aggregate.destination     = sum(perAgent.destination);
    aggregate.totalDifference = sum(perAgent.totalDifference);

    decomposition.perAgent  = perAgent;
    decomposition.aggregate = aggregate;
end

