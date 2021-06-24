classdef oop_stateclass
    
    % state class contains intra-state epidemic SIR model and its solver. 
    
    properties
        eps
        I0
        R0
        Ti
        Tf
        X
        t
        
        nu  % recovery rate
        beta  % infection rate
    end
    
    properties (Dependent)
        S0
        sigma
        discr
        dXdt
        
        % results
        exp_imax
        s_inf
        real_imax
        
        % matrix form
        intra_matrix
    end
    
    
    methods
        % Set S0
        function S0 = get.S0(obj)
            S0 = 1-obj.eps;
        end
        
        % Set sigma
        function sigma = get.sigma(obj)
            sigma = obj.beta/obj.nu;
        end
        
        % Set discriminant
        function discr = get.discr(obj)
            discr = obj.sigma*obj.S0;
        end
        
        % Set dXdt
        function dXdt = get.dXdt(obj)
            function dXdt = eqns(t,X,beta,nu)
                dS = -beta*X(1)*X(2);
                dI = beta*X(1)*X(2)-nu*X(2);
                dR = nu*X(2);
                dXdt = [dS;dI;dR];
            end
            dXdt = @(t,X) eqns(t, X, obj.beta, obj.nu);
        end
        
        % Get the solution X
        function Xval = get.X(obj)
            [~, Xval,~,~,~] = obj.SIR_RKODE45Epi;
        end
        
        % Get the time step
        function tval = get.t(obj)
            [tval, ~,~,~,~] = obj.SIR_RKODE45Epi;
        end
               
        % Solver
        function [t, soln, exp_imax, s_inf, real_imax]=SIR_RKODE45Epi(obj)
            
            if (obj.discr>1)
             a = 0;
             b = 1/obj.sigma;

%              obj.s_inf_plt = plot(x,obj.I0+obj.S0-x+log(x/obj.S0)/obj.sigma,'k','LineWidth',3);
%              grid
             
             exp_imax = obj.I0+obj.S0-1/obj.sigma-(log(obj.sigma*obj.S0))/obj.sigma
             s_inf = fzero(@(x)(obj.I0+obj.S0-x+log(x/obj.S0)/obj.sigma),[a+1.e-9, b])
            end

            % RK_ODE45 step
            options = odeset('RelTol',1e-8,'AbsTol',1e-10);
            [t,soln] = ode45(obj.dXdt,[obj.Ti,obj.Tf],[obj.S0;obj.I0;obj.R0],options);

            real_imax = max(soln(:, 2))

            t_size = size(t);
            soln_size = size(soln)
        end
        
        %%%% Plot things
        function plotting(obj)
            options = odeset('RelTol',1e-8,'AbsTol',1e-10);
            
            if (obj.discr>1)
                a = 0;
                b = 1/obj.sigma;
                h = 1.e-4; % increments
                x = a+0.1*1.e-8:1.e-10:a+0.5*1.e-8;

                figure(1), clf
                plot(x,obj.I0+obj.S0-x+log(x/obj.S0)/obj.sigma,'k','LineWidth',3);
                grid
            end
             
            figure(2),clf
            plot(obj.t,obj.X(:,1),'x-','LineWidth',2)
            hold on 
            plot(obj.t,obj.X(:,2),'x-','LineWidth',2)
            hold on 
            plot(obj.t,obj.X(:,3),'x-','LineWidth',2)
            hold on
            legend({'Susceptibles','Infected','Recovered'},'Location','northeast')

            Int1=cumtrapz(obj.t,obj.X(:,2));
            Int2=trapz(obj.t,obj.X(:,2));

            figure(3),clf
            plot(obj.t, Int1,'LineWidth',2)
            title("Total Infected over time", 'FontSize',18);

            % Plot phase portrait

            % setup
            LW = 2;

            % First contruct a vector of different I0's
            I0_vector = 0:0.05:1;
            S0_vector = 1-I0_vector;

            % Plot S vs. I for each I0
            figure(4), clf, hold on
            axis square;
            axis([0,1,0,1]);
            for i=1:length(I0_vector)
                [tval,Xval] = ode45(obj.dXdt,[obj.Ti,obj.Tf],[S0_vector(i);I0_vector(i);0],options);
                pp= plot(Xval(:,1),Xval(:,2),'LineWidth',LW);
                plot_dir(Xval(1:20:end,1),Xval(1:20:end,2));
            end

            I0_vector = 1.e-6; % if we start at 0, ds/dt = 0 (equil)--> nothing on graph
            S0_vector = 0.5:0.05:0.9;
            for i=1:length(S0_vector)
                [tval,Xval] = ode45(obj.dXdt,[obj.Ti,obj.Tf],[S0_vector(i);I0_vector;0],options);
                pp2= plot(Xval(:,1),Xval(:,2),'LineWidth',LW);
                plot_dir(Xval(1:20:end,1),Xval(1:20:end,2));
            end

            % Plot the line in the diagonal
            xvals = linspace(0,1,200);
            yvals = 1-xvals;
            figure(4), diagline = plot(xvals, yvals,'k','LineWidth',LW);

            % Plot s_max
            s_max = 1/obj.sigma;
            figure(4), maxpt = plot(s_max, 0, 'ro','MarkerSize', 10);

            legend([pp,diagline, maxpt],"discr = "+ obj.discr, "Diagonal", "s_{max} = ^1/_{\sigma}");
            xlabel("Susceptible fraction, s",'FontSize', 15);
            ylabel("Infected fraction, i",'FontSize',15);

            hold off
        end

    end 
    
    
end
