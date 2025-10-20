function [lower_bound, upper_bound] = SetParameterBounds()
%SETPARAMETERBOUNDS Define admissible ranges for the calibratable parameters.
%   The bounds are used to map parameters into the unconstrained space of the
%   Gauss–Newton search. Scalars apply uniformly across destinations, while
%   vector entries specify location-specific limits.

    % Lower bounds guarantee that productivity and policy parameters remain in
    % economically meaningful ranges.
    lower_bound = struct();
    lower_bound.A            = 1.0;
    lower_bound.B            = 0.5;
    lower_bound.tilde_ttau   = [3; 0.1; 1; 1; 3];
    lower_bound.hat_ttau     = 0.5;
    lower_bound.baseUP_eta   = 0.02;
    lower_bound.baseUP_psi   = 0.02;
    lower_bound.ggamma       = 0.5;
    lower_bound.cchi         = 0.05;

    % Upper bounds prevent the search from venturing into implausible regions
    % that would destabilise the simulation.
    upper_bound = struct();
    upper_bound.A            = 3.5;
    upper_bound.B            = 1.6;
    upper_bound.tilde_ttau   = [8; 2; 4; 6; 10];
    upper_bound.hat_ttau     = 10;
    upper_bound.baseUP_eta   = 0.25;
    upper_bound.baseUP_psi   = 0.25;
    upper_bound.ggamma       = 4.0;
    upper_bound.cchi         = 0.35;

end
