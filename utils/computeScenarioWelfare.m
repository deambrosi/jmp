function welfare = computeScenarioWelfare(vf_path, agentData, params, dims, welfareHorizon)
% COMPUTESCENARIOWELFARE Compute discounted welfare measures for a scenario.
%
%   welfare = computeScenarioWelfare(vf_path, agentData, params, dims,
%                                    welfareHorizon) implements the welfare
%   definition provided by the user: realised value functions are evaluated
%   along each simulated trajectory and discounted with β over the first
%   \tilde{T} periods selected in the main script.
%
%   INPUTS:
%       vf_path        - Cell array of value-function structs for t = 1,...,T
%                        containing fields .V (out of network) and .Vn (in
%                        network)
%       agentData      - Struct with trajectories (.location, .wealth, .state,
%                        .network) returned by the simulation routines
%       params         - Parameter struct providing params.bbeta
%       dims           - Dimension struct (used for bounds checking)
%       welfareHorizon - Scalar \tilde{T} specifying the welfare horizon
%
%   OUTPUT:
%       welfare - Struct with fields:
%                   .valuePanel     [Nagents x \tilde{T}] realised V_i(t)
%                   .discountedPanel[ Nagents x \tilde{T}] β^t V_i(t)
%                   .perAgentTotal  [Nagents x 1] discounted sums W_i
%                   .total          Scalar aggregate W
%                   .betaPowers     Row vector of β^t used in the sums
%                   .horizon        The \tilde{T} applied
%
%   AUTHOR: ChatGPT
%   DATE:   May 2025
% =========================================================================

    requiredFields = {'location','wealth','state','network'};
    for k = 1:numel(requiredFields)
        if ~isfield(agentData, requiredFields{k})
            error('agentData is missing required field "%s".', requiredFields{k});
        end
    end

    numAgents = size(agentData.location, 1);
    Tsim      = numel(vf_path);
    Ttilde    = min(welfareHorizon, Tsim);

    valuePanel     = zeros(numAgents, Ttilde);
    discountedPanel= zeros(numAgents, Ttilde);

    betaPowers = params.bbeta .^ (1:Ttilde);

    for t = 1:Ttilde
        vf_t = vf_path{t};
        loc_t = agentData.location(:, t);
        wea_t = agentData.wealth(:, t);
        sta_t = agentData.state(:, t);
        net_t = agentData.network(:, t);

        linInd = sub2ind([dims.S, dims.Na, dims.N], sta_t, wea_t, loc_t);

        values = vf_t.V(linInd);
        values(net_t == 1) = vf_t.Vn(linInd(net_t == 1));

        valuePanel(:, t) = values;
        discountedPanel(:, t) = betaPowers(t) * values;
    end

    perAgentTotal = sum(discountedPanel, 2);

    welfare.valuePanel      = valuePanel;
    welfare.discountedPanel = discountedPanel;
    welfare.perAgentTotal   = perAgentTotal;
    welfare.total           = sum(perAgentTotal);
    welfare.betaPowers      = betaPowers;
    welfare.horizon         = Ttilde;
end

