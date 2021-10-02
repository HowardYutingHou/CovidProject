% Load data
fake_inter = readtable('fake_inter.csv');
fake_data = readtable('fake_data.csv');
% note: cannot perform operations on individual cells directly. Rather, get
% the columns and then get individual element. eg. fake_inter.s_o(3,1)*8

% Start count CPU time
tStart = cputime;

% Setup
N_c = 3;
N_r = 7;

t0 = 0;
tf = 500;
dt = 1;

% First, create an empty vector whose elements are individual states (objs of
% oop_stateclass)
obj_regs(N_r, 1) = oop_stateclass();

% Assign name and SIR (intra-state) parameters to the objects
for i = 1:length(obj_regs)
    obj_regs(i).name = fake_data.name(i);
    obj_regs(i).S = fake_data.S(i);
    obj_regs(i).I = fake_data.I(i);
    obj_regs(i).R = fake_data.R(i);
    obj_regs(i).mu = fake_data.mu(i);
    obj_regs(i).beta = fake_data.beta(i);
    obj_regs(i).nu = fake_data.gamma(i);
    
    % set lockdown start & end day
    obj_regs(i).ld0 = 0;%fake_data.lockdown_start(i);
    obj_regs(i).ldf = 0;%fake_data.lockdown_end(i);
    
    % set vaccination info
    obj_regs(i).rhom = fake_data.rho_0(i);
    obj_regs(i).rhoM = fake_data.rho_K(i);
    obj_regs(i).a = fake_data.rho_r(i);
    obj_regs(i).tV = fake_data.vac_time(i);
end

%%%%% Solve the problem
[obj_regs, soln] = SIR_solver(N_c, N_r, dt, t0, tf, fake_inter, obj_regs);
[A_L, M, Local, u0] = assembleA_Linear(N_c, N_r, fake_inter, obj_regs, tf);

%%%%% Plot the solution
region = 1; % choose which region to plot

tvals = linspace(t0, tf, (tf-t0)/dt);
Svals = soln((region-1)*N_c +1, :);
Ivals = soln((region-1)*N_c +2, :);
Rvals = soln((region-1)*N_c +3, :);

figure(region), clf, hold on
plot(tvals, Svals, 'x-', 'LineWidth', 2);
plot(tvals, Ivals, 'x-', 'LineWidth', 2);
plot(tvals, Rvals, 'x-', 'LineWidth', 2);
legend({'Susceptibles','Infected','Recovered'},'Location','northeast')
title(obj_regs(region).name, 'FontSize',18);
xlabel('Time (days)');
ylabel('Number');

% Print max infected number 
real_imax = max(Ivals);
real_S_e = Svals(1, end);
real_I_e = Ivals(1, end);
fprintf('Real S_e = %s \nReal I_e = %s\n', real_S_e, real_I_e);

% Thm 2.2 from Math Infections P608
if obj_regs(region).sigma > 1
    obj_regs(region).s_inf = 1/obj_regs(region).sigma*obj_regs(region).Tot;
    obj_regs(region).i_inf = obj_regs(region).mu*(obj_regs(region).sigma-1)/obj_regs(region).beta*obj_regs(region).Tot;
    fprintf('sigma = %s \nShould approach endemic equilibrium with S_e = %s and I_e = %s \n', ...
        obj_regs(region).sigma, obj_regs(region).s_inf, obj_regs(region).i_inf);
else
    fprintf('sigma = %s \nShould approach disease-free equilibrium!', obj_regs(region).sigma);
end


% End count CPU time
tEnd = cputime - tStart;
fprintf('\nCPU time: %s \n', tEnd);


% Get mobility matrices in csv files
writematrix(M, 'inter_state.csv');
writematrix(M+A_L, 'total_mobility.csv');



