
classdef oop_stateclass
    
    % state class contains intra-state endemic SIR model and its solver. 
    
    properties
        name
        
        S
        I
        R
%         Tot
        
        mu % reproduction/death rate
        nu  % recovery rate
        beta  % infection rate/number of adequate contacts
        
        ld0 % lockdown start day
        ldf % lockdown end day
                
        % results (for endemic equil)
        s_inf
        i_inf
    end
    
    
    properties (Dependent)
        Tot % total population
        sigma
    end
    
    
    methods
        function obj = setS(obj,val)
            obj.S = val;
        end
        function obj = setI(obj,val)
            obj.I = val;
        end
        function obj = setR(obj,val)
            obj.R = val;
        end
        function a = get.Tot(obj)
            a = obj.S + obj.I + obj.R;
        end
        
        % Set sigma
        function sigma = get.sigma(obj)
            sigma = obj.beta/(obj.nu+obj.mu);
        end
    end
end