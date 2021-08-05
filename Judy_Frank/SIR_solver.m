function [obj_regs, x] = SIR_solver(N_c, N_r, dt, t0, tf, fake_inter, obj_regs)

t = (tf-t0)/dt;
x = zeros(N_r*N_c,t);

[~, M, ~, u0]= assembleA_Linear(N_c, N_r, fake_inter, obj_regs, 1);
x(:,1) = u0;

% A_NL = assembleA_NL(N_c, N_r, obj_regs, 1);
% (A_L+M)*u0+A_NL

for i = 2:t
    % update S,I,R populations each loop
    temp = M*u0+u0;
    a = 1;
    for b = 1:3:length(temp)
        obj_regs(a).S = temp(b);
        obj_regs(a).I = temp(b+1);
        obj_regs(a).R = temp(b+2);
        a = a+1;
    end

    [A_L, M, ~, u0]= assembleA_Linear(N_c, N_r, fake_inter, obj_regs, i);
    A_NL1 = assembleA_NL(N_c, N_r, obj_regs, i);
    aux1 = (A_L+M)*u0 + A_NL1;
    ustar = u0 + dt.*aux1;
    q=0;
    for p = 1:N_r
        pq = p+q;
        obj_regs(p).S = ustar(pq);
        obj_regs(p).I = ustar(pq+1);
        obj_regs(p).R = ustar(pq+2);
        q = q + N_c - 1;
    end
    A_NL2 = assembleA_NL(N_c, N_r, obj_regs, i);
    aux2 = (A_L+M)*ustar + A_NL2;
    u0 = u0 + dt/2 .* (aux1 + aux2);
    x(:,i) = u0;
end

for i = 1:N_r
    display(obj_regs(i).Tot);
end

