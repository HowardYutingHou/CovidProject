function A_L = assembleA_Linear(N_c, N_r, Local, Mobility)
A_L = zeros(N_c*N_r, N_c*N_r);

j = 0;
for i = 1:N_r
    ii = i + j;
    A_L(ii,ii) = Local(1,i); % alpha (beta is for nonlinear part)
    A_L(ii+1,ii+1) = -Local(3,i); % gamma
    A_L(ii+2,ii+1) = Local(3,i); % gamma
    j = j + N_c - 1;
end


% region 1
% A_L(4,1) = Mobility(1,1);
% A_L(5,2) = Mobility(2,1);
% A_L(6,3) = Mobility(3,1);
% A_L(7,1) = Mobility(1,2);
% A_L(8,2) = Mobility(2,2);
% A_L(9,3) = Mobility(3,2);

k = 1;
for i = 1:N_r-1             % i = 1:2
    for j = 1:N_r           % j = 1:3
        if j == k           % Check if it is the main diagonal
            k = k + N_r;    % k = k + 3
        end
        if k > N_r*N_r      % k must be <= 9
            break
        end
        A_L(k,j) = Mobility(j,i);
        k = k + 1;
    end
end

        
% region 2
% A_L(1,4) = Mobility(1,3);
% A_L(2,5) = Mobility(2,3);
% A_L(3,6) = Mobility(3,3);
% A_L(7,4) = Mobility(1,4);
% A_L(8,5) = Mobility(2,4);
% A_L(9,6) = Mobility(3,4);

k = 1;
for i = N_r:2*(N_r-1)       % i = 3:4
    for j = 1:N_r           % j = 1:3
        j_2 = j + N_r;      % j_2 = 4:6
        if j_2 == k
            k = k + N_r;
        end
         if k > N_r*N_r
            break
        end
        A_L(k,j_2) = Mobility(j,i);
        k = k + 1;
    end
end



% region 3
% A_L(1,7) = Mobility(1,5);
% A_L(2,8) = Mobility(2,5);
% A_L(3,9) = Mobility(3,5);
% A_L(4,7) = Mobility(1,6);
% A_L(5,8) = Mobility(2,6);
% A_L(6,9) = Mobility(3,6);

k = 1;
for i = 2*(N_r-1)+1:2*N_r       % i = 5:6
    for j = 1:N_r               % j = 1:3
        j_3 = j + 2*N_r;        % j_3 = 7:9
        if j_3 == k
            k = k + N_r;
        end
        if k > N_r*N_r
            break
        end

        A_L(k,j_3) = Mobility(j,i);
        k = k + 1;
    end
end
