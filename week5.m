clear; clc; close all;
addpath('../week4');   
mpc = loadcase('case39');
Z = My_makeZ(mpc);
[~, Y] = My_makeY(mpc);

% 验证 YZ ≈ I 
I_est = Y * Z;              
I_ref = eye(size(Y));      
Delta = I_est - I_ref;      

err_F = norm(Delta, 'fro');
err_max = max(abs(Delta(:)));
fprintf('[YZ-I] |·|_max = %.3e, ||·||_F = %.3e\n', err_max, err_F);
