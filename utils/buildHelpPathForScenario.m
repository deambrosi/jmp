function G_path = buildHelpPathForScenario(M_path, params, dims, helpMode)
% BUILDHELPPATHFORSCENARIO Construct help distributions conditional on help mode.
%
%   G_path = buildHelpPathForScenario(M_path, params, dims, helpMode) returns
%   the matrix of help probabilities used when simulating a scenario. When the
%   help mode is set to 'none' the function concentrates all probability mass on
%   the zero-help vector, effectively switching off the network externality.
%
%   INPUTS:
%       M_path  - [N x T] path of migrant masses
%       params  - Parameter struct (requires params.ggamma)
%       dims    - Dimension struct (requires dims.H)
%       helpMode- String: 'endogenous' (default) or 'none'
%
%   OUTPUT:
%       G_path  - [H x T] matrix of help probabilities
%
%   AUTHOR: ChatGPT
%   DATE:   May 2025
% =========================================================================

    T = size(M_path, 2);

    if strcmpi(string(helpMode), "none")
        G_path = zeros(dims.H, T);
        G_path(1, :) = 1;
    else
        G_path = computeG(M_path, params.ggamma);
    end
end

