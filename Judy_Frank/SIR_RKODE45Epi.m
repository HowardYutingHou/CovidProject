xfunction [t,X]=SIR_RKODE45Epi()

eps= 0.1; %(1.25)*(10^(-6));


S0=1-eps;
I0=eps;  
R0=0;    
Ti=0;   
Tf=180;
X=zeros(3,1);     

nu=0.3;
beta=0.9;
sigma = beta/nu;

discr=sigma*S0

if (discr>1)
 a = 0;
 b = 1/sigma;
 h = 1.e-4; % increments
 x = a+0.1*1.e-8:1.e-10:a+0.5*1.e-8;
 figure(1),clf
 plot(x,I0+S0-x+log(x/S0)/sigma,'k','LineWidth',3);
 grid
 exp_imax = I0+S0-1/sigma-(log(sigma*S0))/sigma
 s_inf = fzero(@(x)(I0+S0-x+log(x/S0)/sigma),[a+1.e-9, b]) 
end

options = odeset('RelTol',1e-8,'AbsTol',1e-10);

[t,X]=ode45(@(t,X) SIREpi(t,X,beta,nu),[Ti,Tf],[S0;I0;R0],options);


real_imax = max(X(:, 2))

size(t);
size(X);

figure(2),clf
plot(t,X(:,1),'x-','LineWidth',2)
hold on 
plot(t,X(:,2),'x-','LineWidth',2)
hold on 
plot(t,X(:,3),'x-','LineWidth',2)
hold on
legend({'Susceptibles','Infected','Recovered'},'Location','northeast')

Int1=cumtrapz(t,X(:,2));
Int2=trapz(t,X(:,2));

figure(3),clf
plot(t,Int1,'LineWidth',2)
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
    [t,X]=ode45(@(t,X) SIREpi(t,X,beta,nu),[Ti,Tf],[S0_vector(i);I0_vector(i);0],options);
    pp= plot(X(:,1),X(:,2),'LineWidth',LW);
    plot_dir(X(1:20:end,1),X(1:20:end,2));
end

I0_vector = 1.e-6; % if we start at 0, ds/dt = 0 (equil)--> nothing on graph
S0_vector = 0.5:0.05:0.9;
for i=1:length(S0_vector)
    [t,X]=ode45(@(t,X) SIREpi(t,X,beta,nu),[Ti,Tf],[S0_vector(i);I0_vector;0],options);
    pp2= plot(X(:,1),X(:,2),'LineWidth',LW);
    plot_dir(X(1:20:end,1),X(1:20:end,2));
end

% Plot the line in the diagonal
xvals = linspace(0,1,200);
yvals = 1-xvals;
figure(4), diagline = plot(xvals, yvals,'k','LineWidth',LW);

% Plot s_max
s_max = 1/sigma;
figure(4), maxpt = plot(s_max, 0, 'ro','MarkerSize', 10);

legend([pp,diagline, maxpt],"discr = "+discr, "Diagonal", "s_{max} = ^1/_{\sigma}");
xlabel("Susceptible fraction, s",'FontSize', 15);
ylabel("Infected fraction, i",'FontSize',15);



















hold off
end
