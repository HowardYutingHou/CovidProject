function [A_L, M, Local, u0] = assembleA_Linear(N_c, N_r, fake_inter, obj_regs, t)
A_L = zeros(N_c*N_r, N_c*N_r);
Local = [obj_regs.mu; obj_regs.beta; obj_regs.nu];

% Initial compartmental populations for each region in form of a vector
% [S1; I1; R1; S2; I2; R2; ...]
u0 = zeros(N_c*N_r, 1);
q = 0;
for p = 1:N_r
    pq = p+q;
    u0(pq) = obj_regs(p).S;
    u0(pq+1) = obj_regs(p).I;
    u0(pq+2) = obj_regs(p).R;
    obj_regs(p).rho = obj_regs(p).rhom*exp(obj_regs(p).a*(t-obj_regs(p).tV))/...
        (1+obj_regs(p).rhom/obj_regs(p).rhoM*(exp(obj_regs(p).a*(t-obj_regs(p).tV))-1));
    q = q+N_c-1;
end


%%%%% Intra-state mobility (SIR dynamics)
% 
% Endemic
j = 0;
for i = 1:N_r
    ii = i + j;
    A_L(ii,ii) = -obj_regs(i).mu -obj_regs(i).rho;
    A_L(ii+1,ii+1) = -obj_regs(i).nu - obj_regs(i).mu;
    A_L(ii+2,ii) = obj_regs(i).rho;
    A_L(ii+2, ii+1) = obj_regs(i).nu;
    A_L(ii+2,ii+2) = - obj_regs(i).mu;
    j = j + N_c - 1;
end

% Epidemic
% j = 0;
% for i = 1:N_r
%     ii = i + j;
%     A_L(ii,ii) = -obj_regs(i).rho;
%     A_L(ii+1,ii+1) = -obj_regs(i).nu;
%     A_L(ii+2,ii) = obj_regs(i).rho;
%     A_L(ii+2, ii+1) = obj_regs(i).nu;
%     j = j + N_c - 1;
% end

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