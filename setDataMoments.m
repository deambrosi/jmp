function moments = setDataMoments()
%SETDATAMOMENTS Load survey-based empirical targets for the calibration.
%   The returned structure contains the cross-sectional statistics used to
%   discipline the simulated model moments. All vectors follow the ordering of
%   locations 2–6 (i.e. excluding Venezuela).

    % Preallocate the struct that will collect the empirical targets.
    moments = struct();

    % Average income by destination. Location 2 (Colombian border) is
    % normalised to one so the remaining entries can be interpreted as ratios.
    moments.average_income      =   [1;     2.44;   2.61;   2.82;   4.17];

    % Share of migrants residing in each destination 20 periods after the
    % initial shock.
    moments.share_of_migrants   =   [0.10;  0.175;  0.046;  0.138;  0.046];

    % Fraction of migrants that arrived with help for the locations with
    % survey coverage (2, 3, 4, and 6). Location 5 lacks data and is handled
    % separately in the main calibration script.
    moments.came_with_help      =   [0.31;  0.35;   0.37;   0.45];

    % Share of migrants that travelled directly from Venezuela to each
    % destination without intermediate stops.
    moments.came_directly       =   [1;     0.50;   0.10;   0.65;   0.75];

    % Auxiliary unemployment rates (not currently targeted but kept available
    % for potential extensions). The vectors report the average rate over the
    % first two periods and after ten periods for each destination.
    moments.initial_un_rates    =   [0.28;  0.22;   0.20;   0.16;   0.10];
    moments.final_un_rates      =   [0.15;  0.12;   0.04;   0.06;   0.045];

end
