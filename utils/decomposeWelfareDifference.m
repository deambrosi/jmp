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

    counterPanel  = extractNumericPanel(counter, 'valuePanel');
    baselinePanel = extractNumericPanel(baseline, 'valuePanel');

    numAgents = size(counterPanel, 1);

    if ~isempty(varargin)
        Ttilde = varargin{1};
    elseif isfield(counter, 'horizon')
        Ttilde = counter.horizon;
    else
        Ttilde = size(counterPanel, 2);
    end

    Ttilde = ensureScalarHorizon(Ttilde);

    counterPanel  = padPanelToHorizon(counterPanel, Ttilde, 'counterfactual');
    baselinePanel = padPanelToHorizon(baselinePanel, Ttilde, 'baseline');

    betaVec = beta .^ (1:Ttilde);

    perAgent.leaveEarlier   = zeros(numAgents, 1);
    perAgent.pathQuality    = zeros(numAgents, 1);
    perAgent.destination    = zeros(numAgents, 1);
    perAgent.totalDifference= zeros(numAgents, 1);

    for iAgent = 1:numAgents
        jAgent = matchIdx(iAgent);

        vc = counterPanel(iAgent, 1:Ttilde);
        vb = baselinePanel(jAgent, 1:Ttilde);

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

%% ------------------------------------------------------------------------
function panel = extractNumericPanel(welfareStruct, fieldName)
% EXTRACTNUMERICPANEL Safely obtain a numeric matrix from a welfare struct.

    if ~isfield(welfareStruct, fieldName)
        error('Missing field "%s" in welfare struct.', fieldName);
    end

    candidate = welfareStruct.(fieldName);

    if isnumeric(candidate)
        panel = candidate;
        return;
    end

    if iscell(candidate)
        try
            panel = cell2mat(candidate);
            return;
        catch
            error('Unable to convert cell array field "%s" into numeric matrix.', fieldName);
        end
    end

    if isstruct(candidate)
        fn = fieldnames(candidate);
        for k = 1:numel(fn)
            value = candidate.(fn{k});
            if isnumeric(value)
                panel = value;
                return;
            end
        end
        error('Field "%s" is a struct without numeric subfields.', fieldName);
    end

    error('Unsupported data type (%s) encountered in field "%s".', class(candidate), fieldName);
end

%% ------------------------------------------------------------------------
function horizon = ensureScalarHorizon(candidate)
% ENSURESCALARHORIZON Convert a generic horizon descriptor into a scalar integer.

    horizon = [];

    if isnumeric(candidate) && isscalar(candidate)
        horizon = double(candidate);
    elseif iscell(candidate) && numel(candidate) == 1 && isnumeric(candidate{1}) ...
            && isscalar(candidate{1})
        horizon = double(candidate{1});
    elseif isstruct(candidate)
        fn = fieldnames(candidate);
        for k = 1:numel(fn)
            value = candidate.(fn{k});
            if isnumeric(value) && isscalar(value)
                horizon = double(value);
                break;
            end
        end
    end

    if isempty(horizon)
        error('Unable to interpret welfare horizon of type %s as a scalar numeric.', ...
            class(candidate));
    end

    if ~isfinite(horizon) || horizon < 0
        error('Welfare horizon must be a finite, non-negative scalar.');
    end

    rounded = round(horizon);
    if abs(rounded - horizon) > 1e-8
        warning('Welfare horizon %.10g is not an integer; rounding to %d.', horizon, rounded);
    end

    horizon = rounded;
end

%% ------------------------------------------------------------------------
function panel = padPanelToHorizon(panel, targetHorizon, panelName)
% PADPANELTOHORIZON Extend a value panel by repeating the final observation.

    currentHorizon = size(panel, 2);

    if currentHorizon == 0
        error('Value panel for %s scenario is empty; cannot extend horizon.', panelName);
    end

    if currentHorizon >= targetHorizon
        panel = panel(:, 1:targetHorizon);
        return;
    end

    deficit = round(targetHorizon - currentHorizon);
    if deficit <= 0
        panel = panel(:, 1:targetHorizon);
        return;
    end

    lastColumn = panel(:, currentHorizon);
    padding = repmat(lastColumn, 1, deficit);
    panel = [panel, padding];
end

