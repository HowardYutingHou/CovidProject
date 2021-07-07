function x = SIR_solver(N_c, N_r, u0, dt, t0, tf, fake_inter, obj_regs)

t = (tf-t0)/dt;
x = zeros(N_r*N_c,t);
x(:,1) = u0;

[A_L, M, Local, u0]= assembleA_Linear(N_c, N_r, fake_inter, obj_regs);

u2 = u0 + dt.*((A_L+M)*u0 + assembleA_NL(N_c, N_r, Local, u0, obj_regs));
u1 = u0 + dt/2.*((A_L+M)*u0 + assembleA_NL(N_c, N_r, Local, u0, obj_regs) + (A_L+M)*u2 + assembleA_NL(N_c, N_r, Local, u2, obj_regs));
x(:,2) = u1;

for i = 3:t
    u2 = u1 + dt.*((A_L+M)*u1 + assembleA_NL(N_c, N_r, Local, u1, obj_regs));
    u1 = u1 + dt/2.*((A_L+M)*u1 + assembleA_NL(N_c, N_r, Local, u1, obj_regs) + (A_L+M)*u2 + assembleA_NL(N_c, N_r, Local, u2, obj_regs));
    x(:,i) = u1;
end

