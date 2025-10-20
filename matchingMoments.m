%% matchingMoments.m — Calibrate model parameters to match data moments
%{
% matchingMoments.m orchestrates the block-wise calibration of the spatial
% migration model. The script performs the following high-level steps:
%
% 1. Load the baseline configuration together with survey-based target
%    moments and admissible parameter bounds.
% 2. Simulate the baseline economy to benchmark the initial loss across all
%    empirical moments.
% 3. Iterate over parameter blocks (productivity, migration frictions, help
%    parameters) and, for each block, run a Gauss–Newton search that matches
%    the targeted moments while optionally penalising deviations for
%    auxiliary moment groups.
% 4. Track the outer-loop history until convergence is reached or the
%    maximum number of iterations is exhausted.
%
% Compared with the original implementation, this version replaces the
% derivative-free Nelder–Mead search with a finite-difference Gauss–Newton
% scheme that dramatically reduces the number of model simulations needed to
% align the moments. The penalty structure is now configurable on a per-block
% basis, enabling the income block to ignore all non-income moments and
% letting the migration and help blocks softly discipline each other.
%}

clc; clear; close all;
rng(12345);  % Reproducibility for the Monte Carlo simulation




%cd('~/work/sample_proj/code');
%addpath(genpath(pwd));


fprintf('*** Matching model moments ***\n');

%% 1. Load baseline configuration and data
dims        = setDimensionParam();
baseParams  = SetParameters(dims);
dataRaw     = setDataMoments();
[lbStruct, ubStruct] = SetParameterBounds();
lbStruct    = normalizeBoundFields(lbStruct);
ubStruct    = normalizeBoundFields(ubStruct);

locIndices  = (2:dims.N)';
numLoc      = numel(locIndices);

dataMoments = struct();
dataMoments.average_income    = dataRaw.average_income(:);
dataMoments.share_of_migrants = dataRaw.share_of_migrants(:);
dataMoments.came_directly     = dataRaw.came_directly(:);
dataMoments.came_with_help    = NaN(numLoc, 1);

helpLocations   = [2, 3, 4, 6];                  % No survey data for location 5
helpMask        = ismember(locIndices, helpLocations);
dataMoments.came_with_help(helpMask) = dataRaw.came_with_help(:);

momentWeights = struct();
momentWeights.average_income    = ones(numLoc, 1);
momentWeights.share_of_migrants = ones(numLoc, 1);
momentWeights.came_directly     = ones(numLoc, 1);
momentWeights.came_with_help    = ones(numLoc, 1);
momentWeights.came_with_help(~helpMask) = 0;     % Ignore missing locations

%% 2. Simulation options
simOptions = struct();
simOptions.targetPeriod      = 20;               % Target horizon for cross-sectional comparisons
simOptions.incomeWindow      = 20;               % Rolling window for average income
simOptions.overrideSettings  = struct('Nagents', 3000);

%% 3. Initialize calibratable parameters
calibParams = struct();
calibParams.A            = baseParams.A;
calibParams.B            = baseParams.B;
calibParams.tilde_ttau   = baseParams.tilde_ttau;
calibParams.hat_ttau     = baseParams.hat_ttau;
calibParams.baseUp_eta   = baseParams.baseUp_eta;
calibParams.baseUp_psi   = baseParams.baseUp_psi;
calibParams.ggamma       = baseParams.ggamma;
calibParams.cchi         = baseParams.cchi;

%% 4. Baseline evaluation
fprintf('Simulating baseline configuration...\n');
currentMoments   = simulatedMoments(calibParams, simOptions);
momentNames      = {'average_income','share_of_migrants','came_directly','came_with_help'};
allMomentSpecs   = buildAllMomentSpecs(momentNames);
baselineLoss     = computeTotalLoss(currentMoments, dataMoments, momentWeights, allMomentSpecs);

history = struct('params', calibParams, 'moments', currentMoments, 'loss', baselineLoss);

%% 5. Iterative calibration loop
tolParam        = 1e-3;
maxOuterIter    = 8;
penaltyWeight   = 0.05;

incomeIdx   = 2:numLoc;     % Locations 3–6
shareIdx    = 1:numLoc;     % All non-Venezuelan locations
directIdx   = 1:numLoc;
helpIdx     = find(helpMask);

blockStep2 = struct('name','A','index',3:dims.N);
blockStep3 = [
    struct('name','B','index',3:dims.N);
    struct('name','tilde_ttau','index',1:(dims.N-1));
    struct('name','hat_ttau','index',1:dims.N);
    struct('name','baseUp_eta','index',1);
    struct('name','baseUp_psi','index',1)
    ];
blockStep4 = [
    struct('name','ggamma','index',1);
    struct('name','cchi','index',1)
    ];

targetStep2 = struct('name','average_income','indices',incomeIdx);
targetStep3 = [
    struct('name','share_of_migrants','indices',shareIdx);
    struct('name','came_directly','indices',directIdx)
    ];
targetStep4 = struct('name','came_with_help','indices',helpIdx);

penaltyStep2 = [];  % Income calibration: ignore auxiliary penalties entirely
penaltyStep3 = struct('name','came_with_help','indices',helpIdx,'weight',0.15);  % Encourage stability of help flows
penaltyStep4 = [
    struct('name','share_of_migrants','indices',shareIdx,'weight',0.15);      % Keep migration shares close to the reference
    struct('name','came_directly','indices',directIdx,'weight',0.15)          % ... and discipline direct-migration shares
    ];

for iter = 1:maxOuterIter
    fprintf('\n=== Outer iteration %d ===\n', iter);
    prevParamVec = collectParameterVector(calibParams);

    referenceMoments = currentMoments;
    [calibParams, currentMoments, lossIncome] = optimizeBlock( ...
        calibParams, blockStep2, targetStep2, dataMoments, momentWeights, ...
        lbStruct, ubStruct, simOptions, penaltyWeight, referenceMoments, penaltyStep2);
    fprintf('  Productivity block loss: %.6f\n', lossIncome);

    referenceMoments = currentMoments;
    [calibParams, currentMoments, lossMigrants] = optimizeBlock( ...
        calibParams, blockStep3, targetStep3, dataMoments, momentWeights, ...
        lbStruct, ubStruct, simOptions, penaltyWeight, referenceMoments, penaltyStep3);
    fprintf('  Migration frictions block loss: %.6f\n', lossMigrants);

    referenceMoments = currentMoments;
    [calibParams, currentMoments, lossHelp] = optimizeBlock( ...
        calibParams, blockStep4, targetStep4, dataMoments, momentWeights, ...
        lbStruct, ubStruct, simOptions, penaltyWeight, referenceMoments, penaltyStep4);
    fprintf('  Help parameters block loss: %.6f\n', lossHelp);

    totalLoss = computeTotalLoss(currentMoments, dataMoments, momentWeights, allMomentSpecs);
    fprintf('  Total weighted loss: %.6f\n', totalLoss);

    history(end+1) = struct('params', calibParams, 'moments', currentMoments, 'loss', totalLoss); %#ok<SAGROW>

    paramDiff = max(abs(collectParameterVector(calibParams) - prevParamVec));
    fprintf('  Max parameter change: %.3e\n', paramDiff);

    if paramDiff < tolParam
        fprintf('Convergence tolerance reached. Stopping iterations.\n');
        break;
    end
end

fprintf('\nFinal calibrated parameters:\n');
disp(calibParams);

fprintf('\nFinal simulated moments (columns 2–6 correspond to locations 2–6):\n');
disp(currentMoments);

fprintf('Target data moments (aligned with simulated ordering):\n');
disp(dataMoments);

save('estimationresults.mat');


%% ------------------------------------------------------------------------
%% Local helper functions
function vec = collectParameterVector(params)
%COLLECTPARAMETERVECTOR Stack all calibratable parameters into a vector.
    vec = [params.A(:);
           params.B(:);
           params.tilde_ttau(:);
           params.hat_ttau(:);
           params.baseUp_eta;
           params.baseUp_psi;
           params.ggamma;
           params.cchi];
end

function [paramsOut, momentsOut, targetLoss] = optimizeBlock( ...
        paramsIn, blockSpec, targetSpecs, dataMoments, weights, ...
        lowerBounds, upperBounds, simOptions, penaltyWeight, ...
        referenceMoments, penaltySpecs)
%OPTIMIZEBLOCK Run a Gauss–Newton update for a parameter block.
%   The function keeps the non-targeted parameters fixed, constructs
%   finite-difference gradients for the targeted loss, and applies a damped
%   Gauss–Newton step supplemented by optional penalties that keep auxiliary
%   moments close to the reference simulation.

    % Compose human-readable labels so progress can be streamed to the
    % command window while the expensive simulations are running.
    blockNames = arrayfun(@(s) s.name, blockSpec, 'UniformOutput', false);
    blockLabel = strjoin(unique(blockNames, 'stable'), ', ');

    targetNames = arrayfun(@(s) s.name, targetSpecs, 'UniformOutput', false);
    targetLabel = strjoin(unique(targetNames, 'stable'), ', ');
    if isempty(targetLabel)
        targetLabel = '(none)';
    end

    fprintf('  Optimizing block {%s} targeting {%s}...\n', blockLabel, targetLabel);

    % Map the block parameters to the unconstrained "z" space used by the
    % logistic transform so that Gauss–Newton steps can move freely.
    [lbVec, ubVec] = gatherBoundsVector(paramsIn, blockSpec, lowerBounds, upperBounds);
    theta0          = getBlockVector(paramsIn, blockSpec);
    zCurrent        = invertBoundsVector(theta0, lbVec, ubVec);

    % Evaluate the baseline residuals and keep track of the best candidate
    % encountered during the Gauss–Newton search.
    evaluationCounter = 0;
    [currentEval, success] = evaluateCandidate(zCurrent);
    if ~success
        error('Baseline simulation for block optimisation failed.');
    end

    fprintf(['    Initial loss: total %.6f (target %.6f, penalty %.6f) ' ...
             'using step norm 0.000e+00.\n'], ...
            currentEval.loss, currentEval.targetLoss, currentEval.penaltyLoss);

    maxIter        = 10;
    stepTolerance  = 1e-4;
    damping        = 1e-4;
    fdStep         = 1e-4;

    for k = 1:maxIter
        % Build a finite-difference Jacobian of the combined residual vector
        % with respect to the unconstrained parameters.
        J = computeJacobian(@evaluateCandidate, zCurrent, currentEval.residual, fdStep);

        grad = J' * currentEval.residual;
        normalMatrix = J' * J + damping * eye(size(J,2));
        deltaZ = -normalMatrix \ grad;

        if norm(deltaZ) < stepTolerance
            fprintf('    GN iter %d: step norm %.3e below tolerance; stopping.\n', k, norm(deltaZ));
            break;  % No meaningful improvement expected
        end

        % Line-search the Gauss–Newton direction to secure an improvement in
        % the weighted sum of squared residuals.
        stepSize   = 1;
        improved   = false;
        maxBacktrack = 6;
        fprintf(['    GN iter %d: current total %.6f (target %.6f, penalty %.6f); ' ...
                 'attempting step with norm %.3e.\n'], ...
                k, currentEval.loss, currentEval.targetLoss, currentEval.penaltyLoss, norm(deltaZ));
        for attempt = 1:maxBacktrack
            zTrial = zCurrent + stepSize * deltaZ;
            [trialEval, ok] = evaluateCandidate(zTrial);
            if ok && trialEval.loss < currentEval.loss
                zCurrent   = zTrial;
                currentEval = trialEval;
                improved   = true;
                fprintf(['      Accepted step size %.3f -> total %.6f ' ...
                         '(target %.6f, penalty %.6f).\n'], ...
                        stepSize, currentEval.loss, currentEval.targetLoss, currentEval.penaltyLoss);
                break;
            end
            stepSize = stepSize / 2;
            fprintf('      Backtracking step; new tentative size %.3f.\n', stepSize);
        end

        if ~improved
            % The Gauss–Newton step failed to improve the objective. Exit the
            % loop and keep the current parameters.
            fprintf('    GN iter %d: no improving step found; terminating block optimisation.\n', k);
            break;
        end
    end

    % Map the optimiser solution back to the original parameter space and
    % re-simulate to guarantee consistency.
    newValues  = applyBoundsVector(zCurrent, lbVec, ubVec);
    paramsOut  = setBlockVector(paramsIn, blockSpec, newValues);
    momentsOut = currentEval.moments;
    targetLoss = currentEval.targetLoss;

    % Nested helpers ------------------------------------------------------
    function [result, ok] = evaluateCandidate(zVec)
        % Convert the unconstrained vector back into the original parameter
        % space and simulate the model.
        candidateValues = applyBoundsVector(zVec, lbVec, ubVec);
        candidateParams = setBlockVector(paramsIn, blockSpec, candidateValues);
        result = struct();
        evaluationCounter = evaluationCounter + 1;
        try
            candidateMoments = simulatedMoments(candidateParams, simOptions);
        catch
            ok = false;
            result.residual   = Inf;
            result.loss       = Inf;
            result.targetLoss = Inf;
            result.moments    = [];
            fprintf('      [%s] Simulation %d failed (model error).\n', blockLabel, evaluationCounter);
            return;
        end

        [targetResidual, targetLossLocal] = collectMomentResiduals( ...
            candidateMoments, dataMoments, weights, targetSpecs, 1);

        penaltyResidual = [];
        penaltyLossLocal = 0;
        if ~isempty(penaltySpecs)
            [penaltyResidual, penaltyLossLocal] = collectMomentResiduals( ...
                candidateMoments, referenceMoments, weights, penaltySpecs, penaltyWeight);
        end

        result.moments     = candidateMoments;
        result.targetLoss  = targetLossLocal;
        result.penaltyLoss = penaltyLossLocal;
        result.residual    = [targetResidual; penaltyResidual];
        result.loss        = sum(result.residual.^2);
        fprintf(['      [%s] Simulation %d complete: total %.6f ' ...
                 '(target %.6f, penalty %.6f).\n'], ...
                blockLabel, evaluationCounter, result.loss, ...
                result.targetLoss, result.penaltyLoss);
        ok = true;
    end
end

function J = computeJacobian(fun, basePoint, baseResidual, stepSize)
%COMPUTEJACOBIAN Build a finite-difference Jacobian for Gauss–Newton steps.
    numParams = numel(basePoint);
    numResiduals = numel(baseResidual);
    J = zeros(numResiduals, numParams);
    for j = 1:numParams
        direction = zeros(numParams, 1);
        direction(j) = stepSize;
        [evalForward, okForward] = fun(basePoint + direction);
        if ~okForward
            direction(j) = -stepSize;
            [evalBackward, okBackward] = fun(basePoint + direction);
            if ~okBackward
                % If both perturbations fail, leave the column zeroed so the
                % Gauss–Newton step can still proceed using available
                % directions, and note the issue for the command window log.
                fprintf('        Jacobian column %d: both perturbations failed; column left zeroed.\n', j);
                continue;
            else
                J(:, j) = (evalBackward.residual - baseResidual) / (-stepSize);
            end
        else
            J(:, j) = (evalForward.residual - baseResidual) / stepSize;
        end
    end
end

function [lowerVec, upperVec] = gatherBoundsVector(params, blockSpec, lowerStruct, upperStruct)
%GATHERBOUNDSVECTOR Assemble parameter bounds for the active block.
    lowerVec = [];
    upperVec = [];
    for k = 1:numel(blockSpec)
        field     = blockSpec(k).name;
        idx       = blockSpec(k).index;
        lowField  = expandBoundField(lowerStruct, params, field);
        upField   = expandBoundField(upperStruct,  params, field);
        lowerVec  = [lowerVec; lowField(idx(:))];
        upperVec  = [upperVec; upField(idx(:))];
    end
end

function boundVec = expandBoundField(bounds, params, fieldName)
%EXPANDBOUNDFIELD Broadcast scalar bounds to match parameter dimensions.
    if ~isfield(bounds, fieldName)
        error('Bound for parameter "%s" is missing.', fieldName);
    end
    boundValue = bounds.(fieldName);
    paramValue = params.(fieldName);
    if isscalar(boundValue)
        boundVec = repmat(boundValue, numel(paramValue), 1);
    else
        boundVec = boundValue(:);
    end
end

function vec = getBlockVector(params, blockSpec)
%GETBLOCKVECTOR Extract the subset of parameters that belong to the block.
    vec = [];
    for k = 1:numel(blockSpec)
        field = blockSpec(k).name;
        idx   = blockSpec(k).index;
        value = params.(field);
        if isscalar(value)
            vec = [vec; value];
        else
            vec = [vec; value(idx(:))];
        end
    end
end

function paramsOut = setBlockVector(paramsIn, blockSpec, values)
%SETBLOCKVECTOR Inject updated block values back into the parameter struct.
    paramsOut = paramsIn;
    cursor    = 1;
    for k = 1:numel(blockSpec)
        field  = blockSpec(k).name;
        idx    = blockSpec(k).index;
        count  = numel(idx);
        if isscalar(paramsOut.(field))
            paramsOut.(field) = values(cursor);
        else
            temp = paramsOut.(field);
            temp(idx) = reshape(values(cursor:cursor+count-1), size(temp(idx)));
            paramsOut.(field) = temp;
        end
        cursor = cursor + count;
    end
end

function x = applyBoundsVector(z, lower, upper)
%APPLYBOUNDSVECTOR Map unconstrained search variables to bounded parameters.
    x = lower + (upper - lower) ./ (1 + exp(-z));
end

function z = invertBoundsVector(x, lower, upper)
%INVERTBOUNDSVECTOR Map bounded parameters back to the unconstrained space.
    ratio = (x - lower) ./ (upper - lower);
    ratio = min(max(ratio, 1e-6), 1 - 1e-6);
    z     = log(ratio ./ (1 - ratio));
end

function loss = computeMomentLoss(modelMoments, targetMoments, weights, specs, scale)
%COMPUTEMOMENTLOSS Evaluate the weighted least-squares loss for a set of moments.
    if nargin < 5
        scale = 1;
    end
    residuals = collectMomentResiduals(modelMoments, targetMoments, weights, specs, scale);
    loss = sum(residuals.^2);
end

function [residualVec, loss] = collectMomentResiduals(modelMoments, targetMoments, weights, specs, baseScale)
%COLLECTMOMENTRESIDUALS Stack weighted residuals for use in Gauss–Newton updates.
    if nargin < 5
        baseScale = 1;
    end
    residualVec = zeros(0, 1);
    loss        = 0;
    if isempty(specs)
        return;
    end

    if ~isstruct(specs)
        error('Moment specification must be a struct or struct array.');
    end

    for k = 1:numel(specs)
        name = specs(k).name;
        idx  = [];
        if isfield(specs(k), 'indices')
            idx = specs(k).indices;
        end
        [modelVec, targetVec, weightVec] = extractMomentVectors(modelMoments, targetMoments, weights, name, idx);
        mask = ~isnan(modelVec) & ~isnan(targetVec) & (weightVec > 0);
        if ~any(mask)
            continue;
        end

        localScale = baseScale;
        if isfield(specs(k), 'weight') && ~isempty(specs(k).weight)
            localScale = localScale * specs(k).weight;
        end

        diff       = modelVec(mask) - targetVec(mask);
        localW     = weightVec(mask);
        segment    = sqrt(localScale) .* sqrt(localW) .* diff;
        residualVec = [residualVec; segment(:)]; %#ok<AGROW>
        loss        = loss + sum(segment.^2);
    end
end

function [modelVec, targetVec, weightVec] = extractMomentVectors(model, target, weights, name, indices)
%EXTRACTMOMENTVECTORS Retrieve aligned model, target, and weight vectors.
    modelVecFull  = model.(name);
    targetVecFull = target.(name);
    weightVecFull = weights.(name);
    if isempty(indices)
        indices = 1:numel(modelVecFull);
    end
    modelVec  = modelVecFull(indices);
    targetVec = targetVecFull(indices);
    weightVec = weightVecFull(indices);
end

function specs = buildAllMomentSpecs(allNames)
%BUILDALLMOMENTSPECS Create default specifications for every moment.
    specs = repmat(struct('name','', 'indices',[]), numel(allNames), 1);
    for k = 1:numel(allNames)
        specs(k).name    = allNames{k};
        specs(k).indices = [];
    end
end

function totalLoss = computeTotalLoss(modelMoments, dataMoments, weights, specs)
%COMPUTETOTALLOSS Convenience wrapper for the aggregate calibration loss.
    totalLoss = computeMomentLoss(modelMoments, dataMoments, weights, specs, 1);
end

function bounds = normalizeBoundFields(bounds)
%NORMALIZEBOUNDFIELDS Harmonise legacy bound names with the parameter struct.
    renamePairs = {'baseUP_eta','baseUp_eta'; 'baseUP_psi','baseUp_psi'};
    for i = 1:size(renamePairs, 1)
        oldName = renamePairs{i,1};
        newName = renamePairs{i,2};
        if isfield(bounds, oldName)
            bounds.(newName) = bounds.(oldName);
            bounds = rmfield(bounds, oldName);
        end
    end
end
