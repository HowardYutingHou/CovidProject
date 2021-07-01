function [A_L, M, Local] = assembleA_Linear(N_c, N_r, fake_inter, obj_regs)
A_L = zeros(N_c*N_r, N_c*N_r);

% Local = [alpha1, alpha2, alpha3, alpha4, alpha5, alpha6; 
%             beta1, beta2, beta3, beta4, beta5, beta6; 
%             nu1, nu2, nu3, nu4, nu5, nu6]
Local = [obj_regs.alpha; obj_regs.beta; obj_regs.nu];

% Initial compartmental populations for each region in form of a vector
% [S1; I1; R1; S2; I2; R2; ...]
u0 = zeros(N_c*N_r, 1);
for i = 1:N_c:length(u0)
    u0(i) = obj_regs(round((i+N_c-1)/3)).S0;
    u0(i+1) = obj_regs(round((i+N_c-1)/3)).I0;
    u0(i+2) = obj_regs(round((i+N_c-1)/3)).R0;
end

j = 0;
for i = 1:N_r
    ii = i + j;
    A_L(ii,ii) = Local(1,i); % alpha (beta is for nonlinear part)
    A_L(ii+1,ii+1) = -Local(3,i); % gamma
    A_L(ii+2,ii+1) = Local(3,i); % gamma
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

% Use struct to match each pair to mobility parameters
% f_i = struct;
% for i = 1:size(mpara, 1)
%     f_i(i).name = pair(i);
%     f_i(i).data = mpara(i, :);
% end

% create vector for outward mobility for all regions -> become main
% diagonal for mobility matrix
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


% % region 1
% % M(1,1) = -Mobility(1,1)-Mobility(1,2)
% % M(2,2) = -Mobility(2,1)-Mobility(2,2)
% % M(3,3) = -Mobility(3,1)-Mobility(3,2)
% % M(4,1) = Mobility(1,1)
% % M(5,2) = Mobility(2,1)
% % M(6,3) = Mobility(3,1)
% % M(7,1) = Mobility(1,2)
% % M(8,2) = Mobility(2,2)
% % M(9,3) = Mobility(3,2)
% 
% k = 1;
% for i = 1:N_r-1             % i = 1:2
%     for j = 1:N_r           % j = 1:3
%         if j == k           % Check if it is the main diagonal
%             M(k,j) = -Mobility(j,i)-Mobility(j,i+1);
%             M(k+1,j+1) = -Mobility(j+1,i)-Mobility(j+1,i+1);
%             M(k+2,j+2) = -Mobility(j+2,i)-Mobility(j+2,i+1);
%             k = k + N_r;
%         end
%         M(k,j) = Mobility(j,i);
%         k = k + 1;
%     end
% end
% 
%         
% % region 2
% % M(1,4) = Mobility(1,3)
% % M(2,5) = Mobility(2,3)
% % M(3,6) = Mobility(3,3)
% % M(4,4) = -Mobility(1,3)-Mobility(1,4)
% % M(5,5) = -Mobility(2,3)-Mobility(2,4)
% % M(6,6) = -Mobility(3,3)-Mobility(3,4)
% % M(7,4) = Mobility(1,4)
% % M(8,5) = Mobility(2,4)
% % M(9,6) = Mobility(3,4)
% 
% k = 1;
% for i = N_r:2*(N_r-1)       % i = 3:4
%     for j = 1:N_r           % j = 1:3
%         j_2 = j + N_r;      % j_2 = 4:6
%         if j_2 == k
%             M(k,j_2) = -Mobility(j,i)-Mobility(j,i+1);
%             M(k+1,j_2+1) = -Mobility(j+1,i)-Mobility(j+1,i+1);
%             M(k+2,j_2+2) = -Mobility(j+2,i)-Mobility(j+2,i+1);
%             k = k + N_r;
%         end
%         M(k,j_2) = Mobility(j,i);
%         k = k + 1;
%     end
% end
% 
% 
% % region 3
% % M(1,7) = Mobility(1,5)
% % M(2,8) = Mobility(2,5)
% % M(3,9) = Mobility(3,5)
% % M(4,7) = Mobility(1,6)
% % M(5,8) = Mobility(2,6)
% % M(6,9) = Mobility(3,6)
% % M(7,7) = -Mobility(1,5)-Mobility(1,6)
% % M(8,8) = -Mobility(2,5)-Mobility(2,6)
% % M(9,9) = -Mobility(3,5)-Mobility(3,6)
% 
% k = 1;
% last = true;
% for i = 2*(N_r-1)+1:2*N_r       % i = 5:6
%     for j = 1:N_r               % j = 1:3
%         j_3 = j + 2*N_r;        % j_3 = 7:9
%         if last
%             M(j_3,j_3) = -Mobility(j,i)-Mobility(j,i+1);
%             M(j_3+1,j_3+1) = -Mobility(j+1,i)-Mobility(j+1,i+1);
%             M(j_3+2,j_3+2) = -Mobility(j+2,i)-Mobility(j+2,i+1);
%             last = false;
%         end
%         M(k,j_3) = Mobility(j,i);
%         k = k + 1;
%     end
% end
