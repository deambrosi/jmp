function summary = summarizeWelfareComparison(label, welfareWithNetwork, welfareNoNetwork, networkDecomp, baselineTotal, comparisonDecomp)
% SUMMARIZEWELFARECOMPARISON Assemble key welfare statistics for reporting.
%
%   summary = summarizeWelfareComparison(label, welfareWithNetwork,
%       welfareNoNetwork, networkDecomp, baselineTotal, comparisonDecomp)
%   prepares a compact struct containing totals, percentage changes, and the
%   three-way decomposition for both the network effect and the counterfactual
%   versus benchmark comparison.
%
%   INPUTS:
%       label              - Descriptive scenario label
%       welfareWithNetwork - Struct from computeScenarioWelfare (scenario s)
%       welfareNoNetwork   - Struct from computeScenarioWelfare (scenario \tilde{s})
%       networkDecomp      - Output from decomposeWelfareDifference comparing s
%                            with \tilde{s}
%       baselineTotal      - Scalar benchmark welfare (with network). Pass NaN
%                            when summarising the benchmark itself.
%       comparisonDecomp   - (Optional) decomposition comparing scenario s with
%                            the benchmark. Use [] when not applicable.
%
%   OUTPUT:
%       summary - Struct with fields:
%                   .label
%                   .totalWelfare
%                   .totalWithoutNetwork
%                   .networkGain
%                   .networkPct
%                   .networkBreakdown (levels and shares)
%                   .vsBenchmark (only when comparisonDecomp provided)
%
%   AUTHOR: ChatGPT
%   DATE:   May 2025
% =========================================================================

    summary = struct();
    summary.label               = label;
    summary.totalWelfare        = welfareWithNetwork.total;
    summary.totalWithoutNetwork = welfareNoNetwork.total;
    summary.networkGain         = networkDecomp.aggregate.totalDifference;

    if abs(welfareNoNetwork.total) > eps
        summary.networkPct = 100 * summary.networkGain / welfareNoNetwork.total;
    else
        summary.networkPct = NaN;
    end

    summary.networkBreakdown.levels = networkDecomp.aggregate;
    if abs(networkDecomp.aggregate.totalDifference) > eps
        summary.networkBreakdown.shares.leaveEarlier = networkDecomp.aggregate.leaveEarlier / networkDecomp.aggregate.totalDifference;
        summary.networkBreakdown.shares.pathQuality  = networkDecomp.aggregate.pathQuality  / networkDecomp.aggregate.totalDifference;
        summary.networkBreakdown.shares.destination  = networkDecomp.aggregate.destination  / networkDecomp.aggregate.totalDifference;
    else
        summary.networkBreakdown.shares.leaveEarlier = NaN;
        summary.networkBreakdown.shares.pathQuality  = NaN;
        summary.networkBreakdown.shares.destination  = NaN;
    end

    if nargin >= 6 && ~isempty(comparisonDecomp) && ~isnan(baselineTotal)
        compStruct = struct();
        compStruct.levels = comparisonDecomp.aggregate;
        compStruct.absoluteGain = comparisonDecomp.aggregate.totalDifference;
        if abs(baselineTotal) > eps
            compStruct.percentGain = 100 * comparisonDecomp.aggregate.totalDifference / baselineTotal;
        else
            compStruct.percentGain = NaN;
        end

        if abs(comparisonDecomp.aggregate.totalDifference) > eps
            compStruct.shares.leaveEarlier = comparisonDecomp.aggregate.leaveEarlier / comparisonDecomp.aggregate.totalDifference;
            compStruct.shares.pathQuality  = comparisonDecomp.aggregate.pathQuality  / comparisonDecomp.aggregate.totalDifference;
            compStruct.shares.destination  = comparisonDecomp.aggregate.destination  / comparisonDecomp.aggregate.totalDifference;
        else
            compStruct.shares.leaveEarlier = NaN;
            compStruct.shares.pathQuality  = NaN;
            compStruct.shares.destination  = NaN;
        end

        summary.vsBenchmark = compStruct;
    else
        summary.vsBenchmark = [];
    end
end

