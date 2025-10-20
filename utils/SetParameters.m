function params = SetParameters(dims, overrides)
% SETPARAMETERS Initializes all structural parameters used in the model.
%
%   INPUT:
%       dims   - Struct containing dimension settings (N, K, B, etc.)
%
%   OUTPUT:
%       params - Struct of calibrated model parameters:
%           .bbeta     - Discount factor
%           .ssigma    - Elasticity of intertemporal substitution
%           .CONS      - Scaling constant
%           .A         - Productivity levels by location
%           .B         - Amenities by location
%           .ttau      - Base migration costs (with help discounting)
%           .aalpha    - Help discount factor (0 < alpha < 1)
%           .ggamma    - Elasticity of help probability
%           .nnu       - Variance of idiosyncratic location shock
%           .cchi      - Probability of MIN separation outside Venezuela
%           .G0        - Help offer probabilities when M = 0
%           .Pe        - Markov matrix for productivity (eta)
%           .Pb        - Markov matrix for amenities (psi)
%           .P         - Joint Markov matrix for state s = (eta, psi)
%
%   AUTHOR: Agustin Deambrosi
%   LAST REVISED: April 2025
% =========================================================================% SETPARAMETERS Initialize structural parameters with optional overrides.

    if nargin < 2 || isempty(overrides)
        overrides = struct();
    end

    %% Preferences
    params.bbeta  = 0.996315;
    params.ssigma = 2.00;
    params.CONS   = 1e2;

    %% Location-specific features
    params.A = [0.8; 1.2; 1.8; 2.1; 2.2; 3.2];
    params.B = [1.5; 1.0; 1.2; 1.0; 1.0; 0.8];

    %% Migration cost primitives
    params.tilde_ttau = [5; 0.5; 2; 3; 6];
    params.hat_ttau   = [3; 1; 1; 1; 1; 9];
    params.aalpha     = 0.5;
    params.nnu        = 0.1;

    %% Network-related parameters
    params.ggamma     = 0.25;
    params.cchi       = 0.50;

    %% Markov shocks
    params.baseUp_eta = 0.08;
    params.baseUp_psi = 0.05;
    shockProb         = 0.05;
    topDiff           = 1;

    %% Apply overrides when provided
    overrideFields = fieldnames(overrides);
    for i = 1:numel(overrideFields)
        field = overrideFields{i};
        params.(field) = overrides.(field);
    end

    params.A          = params.A(:);
    params.B          = params.B(:);
    params.tilde_ttau = params.tilde_ttau(:);
    params.hat_ttau   = params.hat_ttau(:);

    %% Build migration cost tensor with help discounting
    tauBase = zeros(dims.N);
    for ii = 1:dims.N
        for jj = ii+1:dims.N
            transport = sum(params.tilde_ttau(ii:jj-1));
            tauBase(ii, jj) = params.hat_ttau(jj) + transport;
            tauBase(jj, ii) = params.hat_ttau(ii) + transport;
        end
    end
    params.ttau = getMigrationCostsWithHelp(tauBase, params.aalpha);

    %% Markov transitions
    params.Pe = buildMarkovMatrixDifficultAtTop(dims.K, params.baseUp_eta, shockProb, topDiff);
    params.Pb = buildMarkovMatrixDifficultAtTop(dims.B, params.baseUp_psi, 0,        topDiff);
    params.P  = kron(params.Pe, params.Pb);

    %% Help probabilities when network is empty
    params.G0 = computeG(zeros(dims.N, 1), params.ggamma);
end

