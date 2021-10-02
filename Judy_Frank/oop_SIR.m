classdef oop_SIR
    % A parent class of the basic SIR model. Subclasses include
    % oop_stateclass (individual states) and/or more complicated models
    % such as SEIR, SEIRD...
    
    properties
        name
        
        S
        I
        R
        
        mu % reproduction/death rate
        nu  % recovery rate
        beta  % infection rate/number of adequate contacts
        
        rho % vaccination rate (follows logistic growth)
        rhom
        rhoM
        a % steepness of logistic curve
        tV % vaccination start time
        
        ld0 % lockdown start day
        ldf % lockdown end day
                
        % results (for endemic equil)
        s_inf
        i_inf
    end










end
