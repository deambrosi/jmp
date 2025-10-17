function params = SetParameters(dims)
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
% =========================================================================

    %% Preferences
    params.bbeta	= 0.996315;    % Discount factor (quarterly)
    params.ssigma	= 2.00;        % CRRA utility curvature
    params.CONS		= 1e2;         % Scaling constant for calibration

    %% Location-specific features
    params.A		= [0.8; 1.2; 1.8; 2.1; 2.2; 3.2]; % Productivity by location (wage shifters)
    params.B		= [1.5; 1; 1.2; 1; 1; 0.8]; % Amenity levels by location
   
    

    %% Migration & choice frictions
    % Base migration cost components.
    % Step-by-step transportation costs between consecutive locations (length N-1).
    params.tilde_ttau = [5; 0.5; 2; 3; 6];

    % One-time settlement costs by destination (length N).
    params.hat_ttau   = [3; 1; 1; 1; 1; 9];

    % Construct symmetric migration cost matrix using the proposed structure.
    params.ttau = zeros(dims.N);
    for ii = 1:dims.N
        for ll = ii+1:dims.N
            transportCost = sum(params.tilde_ttau(ii:ll-1));
            params.ttau(ii, ll) = params.hat_ttau(ll) + transportCost;
            params.ttau(ll, ii) = params.hat_ttau(ii) + transportCost;  % enforce symmetry
        end
    end
    
    %params.ttau		=   [0, 6, 15;    % Migration cost matrix (asymmetric)
    %                       6, 0, 10;
    %                       15, 10, 0];
    params.nnu		= 0.1;       % Scale of i.i.d. location taste shocks
    
    %% Help mechanics (network effects)
    params.aalpha	= 0.5;         % Fractional cost with help
    params.ggamma	= 2;           % Elasticity of help probability
    params.ttau		= getMigrationCostsWithHelp(params.ttau, params.aalpha); 
    params.cchi		= 0.25;        % Probability of leaving MIN
    params.G0		= computeG(zeros(dims.N,1), params.ggamma); 
                                 % Help vector is all-zero → no help offers
    
    %% Markov transitions
    shockProb		= 0.05;        % Negative productivity shock probability
    topDiff			= 1;           % Harder to move up at top of the ladder
    baseUp_eta		= 0.08;        % Upward mobility base for productivity
    baseUp_psi		= 0.05;        % Upward mobility base for amenities

    params.Pe		= buildMarkovMatrixDifficultAtTop(dims.K, baseUp_eta, shockProb, topDiff);
    params.Pb		= buildMarkovMatrixDifficultAtTop(dims.B, baseUp_psi, 0, topDiff);
    params.P		= kron(params.Pe, params.Pb);  % Joint transition over s = (eta, psi)

end
