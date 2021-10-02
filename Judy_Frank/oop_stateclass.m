
classdef oop_stateclass < oop_SIR
    
    % state class contains intra-state endemic SIR model and its solver. 
    
%     properties
%         name
%         
%         S
%         I
%         R
%         
%         mu % reproduction/death rate
%         nu  % recovery rate
%         beta  % infection rate/number of adequate contacts
%         
%         rho % vaccination rate (follows logistic growth)
%         rhom
%         rhoM
%         a % steepness of logistic curve
%         tV % vaccination start time
%         
%         ld0 % lockdown start day
%         ldf % lockdown end day
%                 
%         % results (for endemic equil)
%         s_inf
%         i_inf
%     end
    
    
    properties (Dependent)
        Tot % total population
        sigma
    end
    
    
    methods
        function a = get.Tot(obj)
            a = obj.S + obj.I + obj.R;
        end
        
        % Set sigma
        function sigma = get.sigma(obj)
            sigma = obj.beta/(obj.nu+obj.mu);
        end
    end
end