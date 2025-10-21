function dims = setDimensionParam()
% SETDIMENSIONPARAM Initializes model dimension parameters.
%
%   OUTPUT:
%       dims - Struct containing model dimensions:
%           .N  - Number of locations
%           .K  - Number of productivity states
%           .B  - Number of relative amenity states
%           .S  - Number of joint (eta, psi) state combinations
%           .H  - Number of possible help offer vectors (2^N)
%           .Na - Number of coarse asset grid points
%           .na - Number of fine asset grid points
%
%   AUTHOR: Agustin Deambrosi
%   LAST REVISED: April 2025

    dims.N		= 4;                % Number of locations
    dims.K		= 6;                % Productivity states
    dims.B		= 6;                % Amenity states
    dims.S		= dims.K * dims.B;  % Total joint state count (eta × psi)
    dims.H		= 2^dims.N;         % Number of help offer combinations
    dims.Na		= 15;               % Coarse asset grid
    dims.na		= 5000;             % Fine asset grid

end

