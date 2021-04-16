r1 = oop_stateclass();
r1.eps=0.1; r1.I0=0.1; r1.R0=0; r1.Ti=0; r1.Tf=1800; r1.X=zeros(3,1); 
r1.nu=0.005; r1.beta=0.02;

r2 = oop_stateclass();
r2.eps=0.1; r2.I0=0.1; r2.R0=0; r2.Ti=0; r2.Tf=1800; r2.X=zeros(3,1); 
r2.nu=0.005; r2.beta=0.02;

r3 = oop_stateclass();
r3.eps=0.1; r3.I0=0.1; r3.R0=0; r3.Ti=0; r3.Tf=1800; r3.X=zeros(3,1); 
r3.nu=0.005; r3.beta=0.02;
 
% SIR parameters 
alpha1 = 0; alpha2 = 0; alpha3 = 0; 
beta1 = r1.beta; beta2 = r2.beta; beta3 = r3.beta;
nu1 = r1.nu; nu2 = r2.nu; nu3 = r3.nu;

Local = [alpha1, alpha2, alpha3; beta1, beta2, beta3; nu1, nu2, nu3];


Nc = 3;
Nr = 3;

% Mobility parameters
% mu21 = 0.1; mu31 = 0.04; mu12 = 0.07; mu32 = 0.09; mu13 = 0.03; mu23 = 0.009;
% phi21 = 0.02; phi31 = 0.003; phi12 = 0.01; phi32 = 0.008; phi13 = 0.002; phi23 = 0.1;
% psi21 = 0.05; psi31 = 0.01; psi12 = 0.03; psi32 = 0.009; psi13 = 0.03; psi23 = 0.03;

% Randomly generate mobility matrix.
s = rng;
blah = rand(3,6)*0.1;
rng(s);
Mobility = rand(3,6)*0.1;




















