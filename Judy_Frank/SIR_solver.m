function x = SIR_solver(N_c, N_r, Local, u0, dt, t0, tf, fake_inter, obj_regs)

t = (tf-t0)/dt;
x = zeros(N_r*N_c,t);
x(:,1) = u0;

[A_L, M]= assembleA_Linear(N_c, N_r, fake_inter, obj_regs);

u2 = u0 + dt.*(A_L*u0 + M*u0 + assembleA_NL(N_c, N_r, Local, u0));
u1 = u0 + dt/2.*(A_L*u0 + M*u0 + assembleA_NL(N_c, N_r, Local, u0) + A_L*u2 + M*u2 + assembleA_NL(N_c, N_r, Local, u2));
x(:,2) = u1;

for i = 3:t
    u2 = u1 + dt.*(A_L*u1 + M*u1 + assembleA_NL(N_c, N_r, Local, u1));
    u1 = u1 + dt/2.*(A_L*u1 + M*u1 + assembleA_NL(N_c, N_r, Local, u1) + A_L*u2 + M*u2 + assembleA_NL(N_c, N_r, Local, u2));
    x(:,i) = u1;
end

