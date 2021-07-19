function [A_NL] = assembleA_NL(N_c, N_r, Local, Region, obj_regs, t)
%
% This function sets up the non-linear portion of the matrix A for the
% epidemic SIR model.
%
% Variables: N_c - number of compartments;
%            N_r - number of regions
%            Local - matrix that contains all SIR parameters
%            A_NL - vector, non-linear portion of the system, structured as
%                   [-beta1*S1*I1; beta1*S1*I1; 0; 
%                    -beta2*S2*I2; beta2*S2*I2; 0; 
%                         ...    ;      ...   ; 0]
%            Region - vector of SIR [S1;I1;R1;S2;I2;R2;...] which is the
%                     same as u


A_NL = zeros(N_c*N_r, 1);

obj_regs.Tot

for i = 1:N_r
    j = N_c*(i-1)+1;

    % change beta if the region is in lockdown (assume beta is reduced by
    % 90% with lockdown
    A_NL(j) = -Local(2, i)* (1-0.9*(t > obj_regs(i).ld0 & t < obj_regs(i).ldf)) ...
        * Region(j)*Region(j+1)/obj_regs(i).Tot + Local(1, i)*obj_regs(i).Tot; % -beta*I*S/N + mu*N
    A_NL(j+1) = Local(2, i)* (1-0.9*(t > obj_regs(i).ld0 & t < obj_regs(i).ldf))* Region(j)*Region(j+1)/obj_regs(i).Tot; % -beta*I*S/N
    A_NL(j+2) = 0;
end


