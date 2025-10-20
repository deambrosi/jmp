function [lower_bound, upper_bound] = SetParameterBounds()

	lower_bound = struct();
	lower_bound.A =	1; 	% upper bound for all A
	lower_bound.B = 0.5;	% upper bound for all B
	lower_bound.tilde_ttau 	= [3; 0.1; 1; 1; 3]; % upper bound for each 
	lower_bound.hat_ttau 	= 0.5; % upper bound for all 
	lower_bound.baseUP_eta = 0.02;
	lower_bound.baseUP_psi = 0.02;
	lower_bound.ggamma 	   = 0.5;
	lower_bound.cchi 	   = 0.05;
	
	upper_bound = struct();
	upper_bound.A =	3.5; 	% upper bound for all A
	upper_bound.B = 1.6; 	% upper bound for all B
	upper_bound.tilde_ttau 	= [8; 2; 4; 6; 10]; % upper bound for each 
	upper_bound.hat_ttau 	= 10; % upper bound for all 
	upper_bound.baseUP_eta = 0.25;
	upper_bound.baseUP_psi = 0.25;
	upper_bound.ggamma 	   = 4;
	upper_bound.cchi 	   = 0.35;

end