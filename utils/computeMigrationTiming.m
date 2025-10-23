function timing = computeMigrationTiming(agentData, welfareHorizon)
% COMPUTEMIGRATIONTIMING Derive key migration timing statistics for agents.
%
%   timing = computeMigrationTiming(agentData, welfareHorizon) inspects the
%   simulated trajectories stored in agentData and returns the first exit time
%   from Venezuela (location = 1) and the travel time to the final destination
%   used in the welfare decomposition. Periods beyond the analysis horizon are
%   ignored so that the decomposition aligns with the welfare measure
%   constructed in Main.m.
%
%   INPUTS:
%       agentData      - Struct with fields .location capturing agents' paths
%       welfareHorizon - Scalar \tilde{T} specifying the horizon used when
%                        summing realised values
%
%   OUTPUT:
%       timing - Struct with fields:
%                   .tStar        First period t where location ≠ 1
%                   .tDoubleStar  Number of periods required to reach the
%                                  location occupied at \tilde{T}
%                   .everLeft     Logical flag equal to true if the agent
%                                  leaves Venezuela within the horizon
%
%   AUTHOR: ChatGPT
%   DATE:   May 2025
% =========================================================================

    if ~isfield(agentData, 'location')
        error('agentData must contain a field named "location".');
    end

    locations = agentData.location;
    numAgents = size(locations, 1);
    Ttilde    = min(welfareHorizon, size(locations, 2));

    tStar       = (Ttilde + 1) * ones(numAgents, 1);
    tDoubleStar = zeros(numAgents, 1);
    everLeft    = false(numAgents, 1);

    for iAgent = 1:numAgents
        path = locations(iAgent, 1:Ttilde);
        leaveIdx = find(path ~= 1, 1, 'first');

        if ~isempty(leaveIdx)
            tStar(iAgent) = leaveIdx;
            everLeft(iAgent) = true;

            finalLocation = path(end);
            postLeavePath = path(leaveIdx:Ttilde);
            arrivalIdx = find(postLeavePath == finalLocation, 1, 'first');
            if isempty(arrivalIdx)
                arrivalIdx = numel(postLeavePath);
            end
            tDoubleStar(iAgent) = max(0, arrivalIdx - 1);
        else
            % Agent never leaves Venezuela within the analysed horizon.
            tStar(iAgent) = Ttilde + 1;
            tDoubleStar(iAgent) = 0;
            everLeft(iAgent) = false;
        end
    end

    timing.tStar       = tStar;
    timing.tDoubleStar = tDoubleStar;
    timing.everLeft    = everLeft;
    timing.horizon     = Ttilde;
end

