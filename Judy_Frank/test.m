N_c = 3;
N_r = 6;

t0 = 0;
tf = 180;
dt = 1;

% normalized populations
gaS = 30000/30100; gaI = 100/30100; gaR = 0;
flS = 45000/45300; flI = 300/45300; flR = 0;
alS = 78000/78600; alI = 600/78600; alR = 0;
caS = 100000/100400; caI = 400/100400; caR = 0;
iaS = 62000/62500; iaI = 500/62500; iaR = 0;
txS = 120000/120500; txI = 500/120500; txR = 0;


GA = oop_stateclass();
GA.I0=0.1; GA.R0=0; GA.Ti=0; GA.Tf=1800; GA.X=zeros(3,1); 
GA.nu=0.005; GA.beta=0.01;

FL = oop_stateclass();
FL.I0=0.2; FL.R0=0; FL.Ti=0; FL.Tf=1800; FL.X=zeros(3,1); 
FL.nu=0.007; FL.beta=0.02;

AL = oop_stateclass();
AL.I0=0.3; AL.R0=0; AL.Ti=0; AL.Tf=1800; AL.X=zeros(3,1); 
AL.nu=0.009; AL.beta=0.03;

CA = oop_stateclass();
CA.I0=0.1; CA.R0=0; CA.Ti=0; CA.Tf=1800; CA.X=zeros(3,1); 
CA.nu=0.005; CA.beta=0.01;

IA = oop_stateclass();
IA.I0=0.2; IA.R0=0; IA.Ti=0; IA.Tf=1800; IA.X=zeros(3,1); 
IA.nu=0.007; IA.beta=0.02;

TX = oop_stateclass();
TX.I0=0.3; TX.R0=0; TX.Ti=0; TX.Tf=1800; TX.X=zeros(3,1); 
TX.nu=0.009; TX.beta=0.03;

% SIR parameters 
alpha1 = 0; alpha2 = 0; alpha3 = 0; alpha4 = 0; alpha5 = 0; alpha6 = 0;
beta1 = GA.beta; beta2 = FL.beta; beta3 = AL.beta; beta4 = CA.beta; beta5 = IA.beta; beta6 = TX.beta;
nu1 = GA.nu; nu2 = FL.nu; nu3 = AL.nu; nu4 = CA.nu; nu5 = IA.nu; nu6 = TX.nu;

Local = [alpha1, alpha2, alpha3, alpha4, alpha5, alpha6; 
            beta1, beta2, beta3, beta4, beta5, beta6; 
            nu1, nu2, nu3, nu4, nu5, nu6]


% Mobility parameters
% mu21 = 0.1; mu31 = 0.04; mu12 = 0.07; mu32 = 0.09; mu13 = 0.03; mu23 = 0.009;
% phi21 = 0.02; phi31 = 0.003; phi12 = 0.01; phi32 = 0.008; phi13 = 0.002; phi23 = 0.1;
% psi21 = 0.05; psi31 = 0.01; psi12 = 0.03; psi32 = 0.009; psi13 = 0.03; psi23 = 0.03;

% Randomly generate mobility matrix.
 s = rng;
 blah = rand(6,36)*0.1;
 rng(s);
 Mobility = rand(6,36)*0.1;

% Mobility matrix



% Initial compartmental populations for each region in form of a vector
Region = [GA.S0; GA.I0; GA.R0; FL.S0; FL.I0; FL.R0; AL.S0; AL.I0; AL.R0;
            CA.S0; CA.I0; CA.R0; IA.S0; IA.I0; IA.R0; TX.S0; TX.I0; TX.R0;];


[A_L,M] = assembleA_L6sb(N_c, N_r, Local, Mobility);
A_NL = assembleA_NL(N_c, N_r, Local, Region);

% Solve the system
soln = SIR_solver(N_c, N_r, Local, Mobility, Region, dt, t0, tf);

% Plot the solution
region = 3; % choose which region to plot

tvals = linspace(t0, tf, (tf-t0)/dt);

Svals = soln((region-1)*N_c +1, :);

Ivals = soln((region-1)*N_c +2, :);

Rvals = soln((region-1)*N_c +3, :);

figure(1), clf, hold on
plot(tvals, Svals, 'x-', 'LineWidth', 2);
plot(tvals, Ivals, 'x-', 'LineWidth', 2);
plot(tvals, Rvals, 'x-', 'LineWidth', 2);
legend({'Susceptibles','Infected','Recovered'},'Location','northeast')
title("Plot for region " + region, 'FontSize',18);












