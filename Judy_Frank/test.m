% Load data
fake_inter = readtable('fake_inter.csv');
fake_data = readtable('fake_data.csv');
% note: cannot perform operations on individual cells directly. Rather, get
% the columns and then get individual element. eg. fake_inter.s_o(3,1)*8

tStart = cputime;
% Setup
N_c = 3;
N_r = 7;

t0 = 0;
tf = 100;
dt = 1;

% First, create an empty vector whose elements are individual states (objs of
% oop_stateclass)
obj_regs(N_r, 1) = oop_stateclass();


% Assign name and SIR (intra-state) parameters to the objects
for i = 1:length(obj_regs)
    obj_regs(i).name = fake_data.name(i);
    obj_regs(i).Tot = fake_data.N(i);
    obj_regs(i).I0 = fake_data.I(i);
    obj_regs(i).R0 = fake_data.R(i);
    obj_regs(i).alpha = fake_data.alpha(i);
    obj_regs(i).beta = fake_data.beta(i);
    obj_regs(i).nu = fake_data.gamma(i);
end


%%%%% Solve the problem
[A_L,M, Local, u0] = assembleA_Linear(N_c, N_r, fake_inter, obj_regs);
A_NL = assembleA_NL(N_c, N_r, Local, u0, obj_regs);

% Solve the system
soln = SIR_solver(N_c, N_r, u0, dt, t0, tf, fake_inter, obj_regs);

% Plot the solution
region = 2; % choose which region to plot

tvals = linspace(t0, tf, (tf-t0)/dt);

Svals = soln((region-1)*N_c +1, :);

Ivals = soln((region-1)*N_c +2, :);

Rvals = soln((region-1)*N_c +3, :);

figure(1), clf, hold on
plot(tvals, Svals, 'x-', 'LineWidth', 2);
plot(tvals, Ivals, 'x-', 'LineWidth', 2);
plot(tvals, Rvals, 'x-', 'LineWidth', 2);
legend({'Susceptibles','Infected','Recovered'},'Location','northeast')
title(obj_regs(region).name, 'FontSize',18);

% Print max infected number 
% real_imax = max(Ivals)
% if obj_regs(region).sigma > 1
%     s_e = 1/obj_regs(region).sigma
%     i_e = obj_regs(region).sigma
% end

tEnd = cputime - tStart

% Get mobility matrices in csv files
writematrix(M, 'inter_state.csv');
writematrix(M+A_L, 'total_mobility.csv');








