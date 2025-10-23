%% Main Script: Evaluate Benchmark and Policy Counterfactuals with Welfare Accounting
% This script performs the following tasks:
%   1. Initializes model parameters, grids, and functional forms.
%   2. Solves for the steady-state value and policy functions assuming no help (G = G0).
%   3. Defines a catalogue of scenarios (benchmark plus counterfactual aid policies).
%   4. For each scenario, solves the dynamic equilibrium, simulates agent histories,
%      and computes discounted welfare over a user-defined horizon \tilde{T} < T.
%   5. Repeats each scenario with network effects switched off to measure W^{\tilde{s}}.
%   6. Implements the three-way welfare decomposition and produces summary tables.
%
% AUTHOR: Agustin Deambrosi (extended by ChatGPT)
% DATE: May 2025
% VERSION: 3.0
% =========================================================================

clc; clear; close all;

%% 1. Initialize Model Parameters, Grids, and Settings
fprintf('Initializing model parameters and grids...\n');
try
    dims                = setDimensionParam();                            % Dimensions (N, K, B, etc.)
    params              = SetParameters(dims);                            % Structural parameters
    [grids, indexes]    = setGridsAndIndices(dims);                       % Grids and index matrices
    matrices            = constructMatrix(dims, params, grids, indexes);  % Precomputed utility and wealth matrices
    settings            = IterationSettings();                            % Iteration controls and simulation length
    m0                  = createInitialDistribution(dims, settings);      % Initial agent distribution
    fprintf('Initialization completed successfully.\n');
catch ME
    error('Error during initialization: %s', ME.message);
end

tic;

%% 2. Solve No-Help Equilibrium (Steady-State)
fprintf('\nComputing No-Help Value and Policy Functions...\n');
try
    [vf_nh, pol_nh] = noHelpEqm(dims, params, grids, indexes, matrices, settings);  % Solve with G = G0
    fprintf('No-Help Value and Policy Functions found successfully.\n');
catch ME
    error('Error computing No-Help equilibrium: %s', ME.message);
end

%% 3. Initialize Guess for Dynamic Network Agent Distribution
% We assume all agents begin in location 1 and are part of the network.
M_init          = zeros(dims.N, 1);
M_init(1)       = 1;
M0              = repmat(M_init, 1, settings.T);  % [N x T] initial guess for M1

%% 4. Welfare Horizon and Scenario Catalogue
analysisHorizon = min(40, settings.T);   % \tilde{T}: welfare horizon shorter than full simulation

transportMass = 0.90 * ones(dims.N, 1);
transportMass(1:2) = 0;                  % No artificial mass in Venezuela or destination 2

scenarioCatalog(1) = struct('name', 'benchmark', ...
    'label', 'Benchmark', 'type', 'benchmark');
scenarioCatalog(2) = struct('name', 'transport_t1', ...
    'label', 'Transport Aid (start t=1)', 'type', 'transport', ...
    'startPeriod', 1, 'wealthThreshold', 13, 'massIncrease', transportMass, 'budget', 3000);
scenarioCatalog(3) = struct('name', 'transport_t6', ...
    'label', 'Transport Aid (start t=6)', 'type', 'transport', ...
    'startPeriod', 6, 'wealthThreshold', 13, 'massIncrease', transportMass, 'budget', 3000);
scenarioCatalog(4) = struct('name', 'shelter_t1', ...
    'label', 'Shelter Aid (start t=1)', 'type', 'shelter', ...
    'startPeriod', 1, 'wealthThreshold', 8, 'transferAmount', 0.7, ...
    'grantProbability', 0.90, 'budget', 3000);
scenarioCatalog(5) = struct('name', 'shelter_t6', ...
    'label', 'Shelter Aid (start t=6)', 'type', 'shelter', ...
    'startPeriod', 6, 'wealthThreshold', 8, 'transferAmount', 0.7, ...
    'grantProbability', 0.90, 'budget', 3000);

numScenarios = numel(scenarioCatalog);
scenarioResults = repmat(struct('spec', [], 'withNetwork', [], 'withoutNetwork', [], ...
    'networkDecomp', [], 'comparisonDecomp', []), numScenarios, 1);

%% 5. Solve and Evaluate Each Scenario
fprintf('\nEvaluating benchmark and counterfactual scenarios...\n');

for sIdx = 1:numScenarios
    spec = scenarioCatalog(sIdx);
    spec.helpMode = 'endogenous';

    seedOffset = (sIdx - 1) * 100;

    fprintf('\n  -> %s (with network effects)\n', spec.label);
    scenarioWith = runScenario(spec, seedOffset, true, M0, vf_nh, m0, dims, params, ...
        grids, indexes, matrices, settings, analysisHorizon);

    fprintf('  -> %s (network effects disabled)\n', spec.label);
    specNoNet = spec;
    specNoNet.helpMode = 'none';
    scenarioWithout = runScenario(specNoNet, seedOffset, false, M0, vf_nh, m0, dims, ...
        params, grids, indexes, matrices, settings, analysisHorizon);

    identityMatch = (1:settings.Nagents)';
    networkDecomp = decomposeWelfareDifference(scenarioWith.welfare, scenarioWithout.welfare, ...
        params.bbeta, identityMatch, scenarioWith.timing, scenarioWithout.timing, analysisHorizon);

    scenarioResults(sIdx).spec            = spec;
    scenarioResults(sIdx).withNetwork     = scenarioWith;
    scenarioResults(sIdx).withoutNetwork  = scenarioWithout;
    scenarioResults(sIdx).networkDecomp   = networkDecomp;
    scenarioResults(sIdx).comparisonDecomp= [];
end

benchmarkResult = scenarioResults(1);

for sIdx = 2:numScenarios
    matchIdx = matchAgentsForDecomposition(m0, benchmarkResult.withNetwork.timing, ...
        scenarioResults(sIdx).withNetwork.timing, settings, sIdx);
    compDecomp = decomposeWelfareDifference(scenarioResults(sIdx).withNetwork.welfare, ...
        benchmarkResult.withNetwork.welfare, params.bbeta, matchIdx, ...
        scenarioResults(sIdx).withNetwork.timing, benchmarkResult.withNetwork.timing, ...
        analysisHorizon);
    scenarioResults(sIdx).comparisonDecomp = compDecomp;
end

%% 6. Summaries and Persistence
baselineTotal = benchmarkResult.withNetwork.welfare.total;

summaries = repmat(struct(), numScenarios, 1);
for sIdx = 1:numScenarios
    spec = scenarioResults(sIdx).spec;
    compDecomp = scenarioResults(sIdx).comparisonDecomp;
    if isempty(compDecomp)
        baseTotal = NaN;
    else
        baseTotal = baselineTotal;
    end
    summaries(sIdx) = summarizeWelfareComparison(spec.label, ...
        scenarioResults(sIdx).withNetwork.welfare, ...
        scenarioResults(sIdx).withoutNetwork.welfare, ...
        scenarioResults(sIdx).networkDecomp, baseTotal, compDecomp);
end

scenarioNames           = {summaries.label}';
totalWelfare            = [summaries.totalWelfare]';
noNetworkWelfare        = [summaries.totalWithoutNetwork]';
networkGain             = [summaries.networkGain]';
networkPercent          = [summaries.networkPct]';
vsBenchmarkGain         = arrayfun(@(s) extractComparisonMetric(s, 'absoluteGain'), summaries);
vsBenchmarkPercent      = arrayfun(@(s) extractComparisonMetric(s, 'percentGain'), summaries);

welfareTable = table(scenarioNames, totalWelfare, noNetworkWelfare, networkGain, ...
    networkPercent, vsBenchmarkGain, vsBenchmarkPercent, ...
    'VariableNames', {'Scenario','TotalW','NoNetworkW','NetworkGain','NetworkPct', ...
    'VsBenchmarkGain','VsBenchmarkPct'});

resultsDir = fullfile('results', 'welfare');
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end

writetable(welfareTable, fullfile(resultsDir, 'welfare_summary.csv'));
save(fullfile(resultsDir, 'scenario_results.mat'), 'scenarioResults', 'summaries', 'analysisHorizon');

fprintf('\nWelfare summary (values in levels and percentages):\n');
disp(welfareTable);

%% 7. Report total runtime
elapsedTime = toc;
fprintf('\nFull script completed in %.2f seconds.\n', elapsedTime);

%% ------------------------------------------------------------------------
function scenarioOutput = runScenario(spec, seedOffset, shouldPlot, M0, vf_terminal, m0, ...
    dims, params, grids, indexes, matrices, settings, analysisHorizon)
% RUNSCENARIO Solve equilibrium and simulate agents for a single scenario.

    rng(settings.rngSeed + seedOffset, 'twister');

    [pol_eqm, M_eqm, iterations, vf_path] = solveDynamicEquilibrium(M0, vf_terminal, ...
        m0, dims, params, grids, indexes, matrices, settings, spec);

    baseHelp = buildHelpPathForScenario(M_eqm, params, dims, spec.helpMode);

    switch lower(string(spec.type))
        case "transport"
            G_aug = buildTransportAidHelpPath(M_eqm, params, spec, baseHelp);
            [M_total, M_network, agentData, stats] = simulateAgentsTransportAid(m0, pol_eqm, ...
                baseHelp, G_aug, dims, params, grids, settings, spec);
        case "shelter"
            [M_total, M_network, agentData, stats] = simulateAgentsShelterAid(m0, pol_eqm, ...
                baseHelp, dims, params, grids, settings, spec);
        otherwise
            [M_total, M_network, agentData] = simulateAgents(m0, pol_eqm, baseHelp, ...
                dims, params, grids, settings);
            stats = struct();
    end

    if shouldPlot
        plotOutcomeCase(M_total, M_network, agentData, dims, settings, spec.name);
    end

    welfare = computeScenarioWelfare(vf_path, agentData, params, dims, analysisHorizon);
    timing  = computeMigrationTiming(agentData, analysisHorizon);

    scenarioOutput = struct('spec', spec, 'policy', pol_eqm, 'M', M_eqm, ...
        'iterations', iterations, 'vf_path', vf_path, 'agentData', agentData, ...
        'stats', stats, 'welfare', welfare, 'timing', timing, ...
        'M_total', M_total, 'M_network', M_network);
end

%% ------------------------------------------------------------------------
function value = extractComparisonMetric(summaryStruct, fieldName)
% EXTRACTCOMPARISONMETRIC Safely obtain a metric from the benchmark comparison.

    if isempty(summaryStruct.vsBenchmark)
        value = NaN;
    else
        value = summaryStruct.vsBenchmark.(fieldName);
    end
end
