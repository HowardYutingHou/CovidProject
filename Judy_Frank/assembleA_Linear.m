function [A_L, M] = assembleA_Linear(N_c, N_r, Local, Mobility)
A_L = zeros(N_c*N_r, N_c*N_r);
M = zeros(N_c*N_r, N_c*N_r);

j = 0;
for i = 1:N_r
    ii = i + j;
    A_L(ii,ii) = Local(1,i); % alpha (beta is for nonlinear part)
    A_L(ii+1,ii+1) = -Local(3,i); % gamma
    A_L(ii+2,ii+1) = Local(3,i); % gamma
    j = j + N_c - 1;
end


% region 1
% M(1,1) = -Mobility(1,1)-Mobility(1,2)
% M(2,2) = -Mobility(2,1)-Mobility(2,2)
% M(3,3) = -Mobility(3,1)-Mobility(3,2)
% M(4,1) = Mobility(1,1)
% M(5,2) = Mobility(2,1)
% M(6,3) = Mobility(3,1)
% M(7,1) = Mobility(1,2)
% M(8,2) = Mobility(2,2)
% M(9,3) = Mobility(3,2)

k = 1;
for i = 1:N_r-1             % i = 1:2
    for j = 1:N_r           % j = 1:3
        if j == k           % Check if it is the main diagonal
            M(k,j) = -Mobility(j,i)-Mobility(j,i+1);
            M(k+1,j+1) = -Mobility(j+1,i)-Mobility(j+1,i+1);
            M(k+2,j+2) = -Mobility(j+2,i)-Mobility(j+2,i+1);
            k = k + N_r;
        end
        M(k,j) = Mobility(j,i);
        k = k + 1;
    end
end

        
% region 2
% M(1,4) = Mobility(1,3)
% M(2,5) = Mobility(2,3)
% M(3,6) = Mobility(3,3)
% M(4,4) = -Mobility(1,3)-Mobility(1,4)
% M(5,5) = -Mobility(2,3)-Mobility(2,4)
% M(6,6) = -Mobility(3,3)-Mobility(3,4)
% M(7,4) = Mobility(1,4)
% M(8,5) = Mobility(2,4)
% M(9,6) = Mobility(3,4)

k = 1;
for i = N_r:2*(N_r-1)       % i = 3:4
    for j = 1:N_r           % j = 1:3
        j_2 = j + N_r;      % j_2 = 4:6
        if j_2 == k
            M(k,j_2) = -Mobility(j,i)-Mobility(j,i+1);
            M(k+1,j_2+1) = -Mobility(j+1,i)-Mobility(j+1,i+1);
            M(k+2,j_2+2) = -Mobility(j+2,i)-Mobility(j+2,i+1);
            k = k + N_r;
        end
        M(k,j_2) = Mobility(j,i);
        k = k + 1;
    end
end


% region 3
% M(1,7) = Mobility(1,5)
% M(2,8) = Mobility(2,5)
% M(3,9) = Mobility(3,5)
% M(4,7) = Mobility(1,6)
% M(5,8) = Mobility(2,6)
% M(6,9) = Mobility(3,6)
% M(7,7) = -Mobility(1,5)-Mobility(1,6)
% M(8,8) = -Mobility(2,5)-Mobility(2,6)
% M(9,9) = -Mobility(3,5)-Mobility(3,6)

k = 1;
last = true;
for i = 2*(N_r-1)+1:2*N_r       % i = 5:6
    for j = 1:N_r               % j = 1:3
        j_3 = j + 2*N_r;        % j_3 = 7:9
        if last
            M(j_3,j_3) = -Mobility(j,i)-Mobility(j,i+1);
            M(j_3+1,j_3+1) = -Mobility(j+1,i)-Mobility(j+1,i+1);
            M(j_3+2,j_3+2) = -Mobility(j+2,i)-Mobility(j+2,i+1);
            last = false;
        end
        M(k,j_3) = Mobility(j,i);
        k = k + 1;
    end
end
