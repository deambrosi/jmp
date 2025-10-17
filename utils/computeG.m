function G = computeG(M, ggamma)
% COMPUTEG Computes the PMF over help offer vectors h ∈ {0,1}^N.
%
%   Given migrant network masses M and elasticity ggamma, this function 
%   computes the joint distribution G(h | M_t) for each possible help vector h, 
%   assuming independent Bernoulli offers with probabilities π^j(M^j_t).
%
%   INPUTS:
%       M       - [N x T] matrix of migrant shares per location over time
%       ggamma  - Elasticity of help probability with respect to migrant mass
%
%   OUTPUT:
%       G       - [2^N x T] matrix where G(h_idx, t) = prob(h = h_idx | M_t)
%
%   AUTHOR: Agustin Deambrosi
%   LAST REVISED: April 2025
% =========================================================================

    [N, T]			= size(M);
    P				= M.^ggamma;                          % π(M) for each t
    Hmat			= dec2bin(0:(2^N - 1)) - '0';         % [2^N x N] all binary help vectors
    G				= zeros(2^N, T);                      % Output: PMF matrix

    for t = 1:T
        pi_t		= P(:, t)';                           % π^j(M^j_t) row vector
        one_minus	= 1 - pi_t;

        for h_idx = 1:2^N
            h		= Hmat(h_idx, :);                    % help vector for this index
            G(h_idx, t) = prod(pi_t.^h .* one_minus.^(1 - h));  % joint pmf value
        end
    end
