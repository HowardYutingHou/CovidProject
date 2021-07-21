function x = SIR_solver(N_c, N_r, u0, dt, t0, tf, fake_inter, obj_regs)

t = (tf-t0)/dt;
x = zeros(N_r*N_c,t);
x(:,1) = u0;

[A_L, M, ~, u0]= assembleA_Linear(N_c, N_r, fake_inter, obj_regs);

% Compute interstate movement and update total populations
temp = M*u0 +u0;
a = 1;
for b = 1:3:length(temp)
    obj_regs(a).setS(0);
%     obj_regs(a).S = temp(b);
    obj_regs(a).I = temp(b+1);
    obj_regs(a).R = temp(b+2);
%     obj_regs(a).Tot = obj_regs(a).S + obj_regs(a).I + obj_regs(a).R;
    a = a+1;
end

% Heun
A_NL = assembleA_NL(N_c, N_r, obj_regs, 1);
u2 = u0 + dt.*((A_L+M)*u0 + A_NL);
u1 = u0 + dt/2.*((A_L+M)*u0 + A_NL + (A_L+M)*u2 + A_NL);
x(:,2) = u1;

for i = 3:t
    % update S,I,R populations each loop
    temp = M*u1+u1;
    a = 1;
    for b = 1:3:length(temp)
        obj_regs(a).S = temp(b);
        obj_regs(a).I = temp(b+1);
        obj_regs(a).R = temp(b+2);
%         obj_regs(a).Tot = obj_regs(a).S + obj_regs(a).I + obj_regs(a).R;
        obj_regs(a).Tot
        a = a+1;
    end
    
    % Heun
    A_NL = assembleA_NL(N_c, N_r, obj_regs, 1);
    u2 = u1 + dt.*((A_L+M)*u1 + A_NL);
    u1 = u1 + dt/2.*((A_L+M)*u1 + A_NL + (A_L+M)*u2 + A_NL);
    x(:,i) = u1;
    
end

for i = 1:N_r
    display(obj_regs(i).Tot);
end