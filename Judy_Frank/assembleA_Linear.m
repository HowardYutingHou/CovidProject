function [A_L, M, Local, u0] = assembleA_Linear(N_c, N_r, fake_inter, obj_regs)
A_L = zeros(N_c*N_r, N_c*N_r);

% Local = [alpha1, alpha2, alpha3, alpha4, alpha5, alpha6; 
%             beta1, beta2, beta3, beta4, beta5, beta6; 
%             nu1, nu2, nu3, nu4, nu5, nu6]
Local = [obj_regs.mu; obj_regs.beta; obj_regs.nu];

% Initial compartmental populations for each region in form of a vector
% [S1; I1; R1; S2; I2; R2; ...]
u0 = zeros(N_c*N_r, 1);
for i = 1:N_c:length(u0)
    u0(i) = obj_regs(round((i+N_c-1)/3)).S0;
    u0(i+1) = obj_regs(round((i+N_c-1)/3)).I0;
    u0(i+2) = obj_regs(round((i+N_c-1)/3)).R0;
end


%%%%% Intra-state mobility (SIR dynamics)
j = 0;
for i = 1:N_r
    ii = i + j;
    A_L(ii,ii) = -Local(1,i); % -mu*S (death rate) (beta is for nonlinear part)
    A_L(ii+1,ii+1) = -Local(3,i) - Local(1,i); % -gamma*I - mu*I
    A_L(ii+2,ii+1) = Local(3,i) - Local(1,i); % gamma - mu*R
    j = j + N_c - 1;
end


%%%%% Mobility matrix

% Concatenate State1 and State2 columns in fake_inter to create pairs for
% keys. Turn 3-8th columns into array to get mobility info in rows.
pair = strcat(fake_inter.State1, {' '}, fake_inter.State2);
pair2 = strcat(fake_inter.State2, {' '}, fake_inter.State1);
mpara = table2array(fake_inter(:, 3:8));
% mpara is in the form [s_i s_o i_i i_o r_i r_o]
% s_i is from State2 to State1; s_o is from State1 to State2

%%% OUTWARD flow (main diagonal)

out1 = zeros(N_c*N_r, 1);
ct = 1;
for st1 = 1:length(obj_regs)-1
    for st2 = 2:length(obj_regs)
        for j = 1:length(pair)
            if strcmpi(strcat(obj_regs(st1).name, {' '}, obj_regs(st2).name), pair(j))
                 out1(ct) = out1(ct) - mpara(j, 2); % use outward mobility (s_o)
                 out1(ct+1) = out1(ct+1) - mpara(j, 4); % i_o
                 out1(ct+2) = out1(ct+2) - mpara(j, 6); % r_o
            end
        end
    end
    ct = ct+3;
end

out2 = zeros(N_c*N_r, 1);
ct = 1+N_c;
for st2 = 2:length(obj_regs)
    for st1 = 1:length(obj_regs)-1
        for j = 1:length(pair)
            if strcmpi(strcat(obj_regs(st2).name, {' '}, obj_regs(st1).name), pair2(j))
                out2(ct) = out2(ct) - mpara(j, 1); % use outward mobility (s_i)
                out2(ct+1) = out2(ct+1) - mpara(j, 3); % i_i
                out2(ct+2) = out2(ct+2) - mpara(j, 5); % r_i
                % disp(pair2(j))
            end
        end
    end
    ct = ct+3;
end

out = out1 + out2; %%%%%% Vector of main diagonal for mobility matrix!!
M = diag(out);

%%% INWARD flow

% First construct portion below main diagonal. Use pair and s_o, i_o, r_o
% info. Fill column-wise.
inw = zeros(N_c*N_r, N_c*N_r);

% individual small diagonal matrix on the x-z plane
pair3d_out = zeros(N_c, length(pair), N_c);
for i = 1: length(pair)
    pair3d_out(:, i, :) = diag([mpara(i, 2); mpara(i, 4); mpara(i, 6)]);
end
pair3d_in = zeros(N_c, length(pair), N_c);
for i = 1: length(pair)
    pair3d_in(:, i, :) = diag([mpara(i, 1); mpara(i, 3); mpara(i, 5)]);
end

ct = 1;
% Below diagonal (fill in column-wise, use pair and s_o, i_o, r_o)
for c = 1:N_c:length(inw)-N_c
    for r = c+N_c:N_c:length(inw)-N_c+1
        if strcmpi(strcat(obj_regs(round((c+N_c)/N_c)).name, {' '}, obj_regs(round((r+N_c)/N_c)).name), pair(ct))
            inw(r:r+2, c:c+2) = pair3d_out(:, ct,:);
        else
            error('Incorrect inter-state movement info!');
        end
        ct = ct+1;
    end
end

ct = 1;
% Above diagonal (fill in row-wise, use pair2 and s_i, r_i, r_I)
for r = 1: N_c: length(inw)-2*N_c+1
    for c = r+N_c:N_c:length(inw)-N_c+1
        if strcmpi(strcat(obj_regs(round((c+N_c)/N_c)).name, {' '}, obj_regs(round((r+N_c)/N_c)).name), pair2(ct))
            inw(r:r+2, c:c+2) = pair3d_in(:, ct,:);
        else
            msg = 'Incorrect inter-state movement info! ' + string(obj_regs(round((c+N_c)/N_c)).name) + ' ' + string(obj_regs(round((r+N_c)/N_c)).name)+ ' ' + string(pair2(ct)); 
            error(msg);
        end
        ct = ct+1;
    end
end

% M = M+inw;

% Test theoremm without mobility
M = zeros(N_c*N_r, N_c*N_r);

