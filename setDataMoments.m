function moments = setDataMoments()


    % create 
    moments = struct();
    
    % average income for all (dims.N - 1) locations that are not Venezuela.
    % locations 2 to 6
    % Colombian border (location 2) is normalized to 1
    moments.average_income      =   [1;     2.44;   2.61;   2.82;   4.17]; 
    
    % share of migrants in each location (exluding venezuela) after 20 periods
    moments.share_of_migrants   =   [0.10;  0.175;  0.046;  0.138;  0.046];
    
    % amount of people that arrive with help (for locations 2, 3, 4, and 6; no data for peru (location 5)):
    moments.came_with_help      =   [0.31;  0.35;   0.37;   0.45];
    
    % amount of people that came directly from location 1 to each destination (without 'settling' in any other location)
    % for all N-1 locations:
    moments.came_directly       =   [1;     0.50;   0.10;   0.65;   0.75];
    
    % unemployment rates for for all (dims.N - 1) locations that are not Venezuela.
    % average unemployment rates at first two periods vs. after 10 periods in the location
    moments.initial_un_rates    =   [0.28;  0.22;   0.20;   0.16;   0.10];
    moments.final_un_rates      =   [0.15;  0.12;   0.04;   0.06;   0.045];

end


