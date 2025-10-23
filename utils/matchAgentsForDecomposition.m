function matchIdx = matchAgentsForDecomposition(m0, baselineTiming, counterTiming, settings, seedOffset)
% MATCHAGENTSFORDECOMPOSITION Construct the rematching map f(i) described by the user.
%
%   matchIdx = matchAgentsForDecomposition(m0, baselineTiming, counterTiming,
%   settings, seedOffset) returns a permutation of agent indices so that agent
%   i in the counterfactual scenario is compared with agent matchIdx(i) in the
%   benchmark scenario when computing welfare differences. The procedure first
%   pairs agents that never leave Venezuela in both scenarios and then matches
%   the remaining agents at random within each initial state group.
%
%   INPUTS:
%       m0             - Initial agent distribution (struct array with fields
%                        .state and .wealth)
%       baselineTiming - Struct from computeMigrationTiming for scenario b
%       counterTiming  - Struct from computeMigrationTiming for scenario c
%       settings       - Iteration settings providing settings.rngSeed
%       seedOffset     - Optional integer added to the base seed to ensure
%                        different counterfactual comparisons produce distinct
%                        yet reproducible matchings
%
%   OUTPUT:
%       matchIdx - [Nagents x 1] vector where matchIdx(i) is the index of the
%                  benchmark agent paired with counterfactual agent i
%
%   AUTHOR: ChatGPT
%   DATE:   May 2025
% =========================================================================

    if nargin < 5 || isempty(seedOffset)
        seedOffset = 0;
    end

    if ~isfield(settings, 'rngSeed')
        error('settings must include an rngSeed field.');
    end

    rng(settings.rngSeed + seedOffset, 'twister');

    numAgents = numel(m0);
    matchIdx  = zeros(numAgents, 1);

    baseKeys = [[m0.state]' [m0.wealth]'];
    [~, ~, groupIds] = unique(baseKeys, 'rows');

    for g = 1:max(groupIds)
        groupMembers = find(groupIds == g);
        basePool     = groupMembers;
        counterPool  = groupMembers;

        baseNever    = basePool(~baselineTiming.everLeft(basePool));
        counterNever = counterPool(~counterTiming.everLeft(counterPool));

        numNeverPairs = min(numel(baseNever), numel(counterNever));
        if numNeverPairs > 0
            baseMatch    = baseNever(1:numNeverPairs);
            counterMatch = counterNever(1:numNeverPairs);
            matchIdx(counterMatch) = baseMatch;

            basePool    = setdiff(basePool, baseMatch, 'stable');
            counterPool = setdiff(counterPool, counterMatch, 'stable');
        end

        if numel(basePool) ~= numel(counterPool)
            error('Mismatch in group sizes after matching non-migrants.');
        end

        permutedBase = basePool(randperm(numel(basePool)));
        matchIdx(counterPool) = permutedBase;
    end
end

