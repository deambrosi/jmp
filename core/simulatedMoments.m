function modelMoments = simulatedMoments(paramOverrides, options)
% SIMULATEDMOMENTS Solve the dynamic equilibrium and extract model moments.
%
%   modelMoments = simulatedMoments(paramOverrides, options)
%
%   paramOverrides : struct with fields to overwrite defaults in SetParameters.
%   options        : struct with optional fields
%                       .targetPeriod       → horizon used for cross-sectional moments
%                       .incomeWindow       → number of periods for average income
%                       .overrideSettings   → struct of IterationSettings overrides
%
%   The function mirrors the baseline pipeline: initialize parameters,
%   solve for the no-help value functions, solve the transition equilibrium,
%   and run a detailed simulation that records the statistics needed for
%   moment matching.

    if nargin < 2 || isempty(options)
        options = struct();
    end

    dims    = setDimensionParam();
    params  = SetParameters(dims, paramOverrides);
    [grids, indexes] = setGridsAndIndices(dims);
    matrices         = constructMatrix(dims, params, grids, indexes);
    settings         = IterationSettings();

    if isfield(options, 'overrideSettings')
        overrideFields = fieldnames(options.overrideSettings);
        for k = 1:numel(overrideFields)
            settings.(overrideFields{k}) = options.overrideSettings.(overrideFields{k});
        end
    end

    m0 = createInitialDistribution(dims, settings);

    [vf_nh, ~] = noHelpEqm(dims, params, grids, indexes, matrices, settings);

    M_init         = zeros(dims.N, 1);
    M_init(1)      = 1;
    M0             = repmat(M_init, 1, settings.T);

    [pol_eqm, M_eqm, ~, ~] = solveDynamicEquilibrium(M0, vf_nh, m0, ...
        dims, params, grids, indexes, matrices, settings, []);

    G_dist = zeros(dims.H, settings.T);
    for t = 1:settings.T
        G_dist(:, t) = computeG(M_eqm(:, t), params.ggamma);
    end

    [~, ~, agentData] = detailedSimulateAgents(m0, pol_eqm, G_dist, dims, params, grids, settings);

    simOpts.targetPeriod = getOption(options, 'targetPeriod', 20);
    simOpts.incomeWindow = getOption(options, 'incomeWindow', simOpts.targetPeriod);

    modelMoments = computeSimulatedMoments(agentData, dims, settings, simOpts);
end

%% ------------------------------------------------------------------------
function value = getOption(options, name, defaultValue)
    if isfield(options, name) && ~isempty(options.(name))
        value = options.(name);
    else
        value = defaultValue;
    end
end

function modelMoments = computeSimulatedMoments(agentData, dims, settings, options)
    targetPeriod = min(options.targetPeriod, settings.T);
    incomeWindow = min(options.incomeWindow, targetPeriod);

    locIndices = (2:dims.N)';
    numLoc     = numel(locIndices);

    location   = agentData.location;
    income     = agentData.income;
    helpArrival= agentData.helpArrival;

    averageIncome = NaN(numLoc, 1);
    for k = 1:numLoc
        loc  = locIndices(k);
        mask = (location(:, 1:incomeWindow) == loc);
        if any(mask(:))
            averageIncome(k) = sum(income(mask)) / sum(mask(:));
        end
    end

    finalLoc      = location(:, targetPeriod);
    migrantMask   = finalLoc ~= 1;
    counts        = accumarray(finalLoc(migrantMask), 1, [dims.N, 1]);
    totalMigrants = sum(counts(2:end));
    shareMigrants = NaN(numLoc, 1);
    if totalMigrants > 0
        shareMigrants = counts(2:end) / totalMigrants;
    end

    cameDirect   = NaN(numLoc, 1);
    cameWithHelp = NaN(numLoc, 1);

    for k = 1:numLoc
        loc       = locIndices(k);
        agentList = find(finalLoc == loc);
        if isempty(agentList)
            continue;
        end

        directCount = 0;
        helpCount   = 0;

        for idx = 1:numel(agentList)
            agentId = agentList(idx);
            path    = location(agentId, 1:targetPeriod);

            firstNonOne = find(path ~= 1, 1, 'first');
            if ~isempty(firstNonOne) && path(firstNonOne) == loc
                directCount = directCount + 1;
            end

            arrivals = find(path(2:targetPeriod) == loc & path(1:targetPeriod-1) ~= loc);
            if ~isempty(arrivals)
                arrivalPeriod = arrivals(end) + 1;
                helpCount     = helpCount + (helpArrival(agentId, arrivalPeriod) > 0);
            end
        end

        cameDirect(k)   = directCount / numel(agentList);
        cameWithHelp(k) = helpCount   / numel(agentList);
    end

    modelMoments.average_income    = averageIncome;
    modelMoments.share_of_migrants = shareMigrants;
    modelMoments.came_directly     = cameDirect;
    modelMoments.came_with_help    = cameWithHelp;
    modelMoments.targetPeriod      = targetPeriod;
    modelMoments.incomeWindow      = incomeWindow;
end
