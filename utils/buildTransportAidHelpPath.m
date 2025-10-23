function G_aug = buildTransportAidHelpPath(M_path, params, scenario, baseHelpPath)
% BUILDTRANSPORTAIDHELPPATH Augment help distributions with transport aid mass.
%
%   G_aug = buildTransportAidHelpPath(M_path, params, scenario, baseHelpPath)
%   replicates the logic used in the counterfactual transport programme. The
%   augmented path coincides with the benchmark path before the programme
%   starts, and—unless helpMode is 'none'—uses the mass-increase vector to shift
%   help probabilities afterwards.
%
%   INPUTS:
%       M_path       - [N x T] equilibrium network masses
%       params       - Parameter struct (uses params.ggamma)
%       scenario     - Scenario struct containing fields massIncrease,
%                      startPeriod, and helpMode
%       baseHelpPath - [H x T] baseline help probabilities
%
%   OUTPUT:
%       G_aug        - [H x T] augmented help probabilities
%
%   AUTHOR: ChatGPT
%   DATE:   May 2025
% =========================================================================

    if ~isfield(scenario, 'massIncrease') || isempty(scenario.massIncrease)
        error('Transportation aid scenario must include a massIncrease field.');
    end

    if ~isfield(scenario, 'startPeriod') || isempty(scenario.startPeriod)
        error('Transportation aid scenario must include a startPeriod field.');
    end

    if nargin < 4 || isempty(baseHelpPath)
        baseHelpPath = computeG(M_path, params.ggamma);
    end

    T = size(M_path, 2);
    startPeriod = max(1, min(T, scenario.startPeriod));

    if strcmpi(string(scenario.helpMode), "none")
        G_aug = baseHelpPath;
        return;
    end

    massIncrease = scenario.massIncrease(:);
    if numel(massIncrease) == 1
        massIncrease = repmat(massIncrease, size(M_path, 1), 1);
    elseif numel(massIncrease) ~= size(M_path, 1)
        error('massIncrease must be a scalar or an N x 1 vector.');
    end

    G_aug = baseHelpPath;
    for t = startPeriod:T
        augmentedMass = M_path(:, t) + massIncrease;
        G_aug(:, t) = computeG(augmentedMass, params.ggamma);
    end
end

