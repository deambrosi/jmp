function shifted = shiftSeriesWithPadding(series, shiftSteps)
% SHIFTSERIESWITHPADDING Shift a time series while padding out-of-range values.
%
%   shifted = shiftSeriesWithPadding(series, shiftSteps) shifts the input
%   SERIES forward (positive shiftSteps) or backward (negative shiftSteps).
%   Values that would fall outside the observed support are padded using the
%   first observation for t ≤ 0 and the last observation for t ≥ T, matching
%   the welfare decomposition convention described in the project notes.
%
%   INPUTS:
%       series     - Row or column vector of length T with values V(t)
%       shiftSteps - Integer number of periods to shift. Positive values move
%                    the series backward in time (evaluate V(t + Δ)), while
%                    negative values move it forward.
%
%   OUTPUT:
%       shifted    - Vector of the same size as SERIES containing the shifted
%                    values with appropriate padding.
%
%   AUTHOR: ChatGPT
%   DATE:   May 2025
% =========================================================================

    if ~isvector(series)
        error('Input series must be a one-dimensional vector.');
    end

    T = numel(series);
    shifted = zeros(size(series));

    for t = 1:T
        sourceIdx = t + shiftSteps;
        if sourceIdx < 1
            shifted(t) = series(1);
        elseif sourceIdx > T
            shifted(t) = series(end);
        else
            shifted(t) = series(sourceIdx);
        end
    end

    % Preserve original orientation (row vs column)
    if iscolumn(series)
        shifted = shifted(:);
    else
        shifted = shifted(:).';
    end
end

