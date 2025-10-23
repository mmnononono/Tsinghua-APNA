clear; clc; close all;
addpath('../week4');   
mpc = loadcase('case39');

% === 计算 Z ===
Z = My_makeZ(mpc);

% === 生成 Ybus（用你自己的 My_makeY） ===
[~, Y] = My_makeY(mpc);

% === 验证 Y^{-1} ≈ Z ===
I  = speye(size(Y,1));
Yinv = full(Y \ I);          % 真正的 Y^{-1}
E = Yinv - Z;                   % 误差矩阵
[err_max, idx] = max(abs(E(:)));
[i_max, j_max] = ind2sub(size(E), idx);

fprintf('\n[Verify] Y^{-1} vs Z：\n');
fprintf('  max|Y^{-1}-Z| = %.3e  @ (i=%d, j=%d)\n', err_max, i_max, j_max);
fprintf('  range |E|     : [%.3e, %.3e]\n', min(abs(E(:))), max(abs(E(:))));

% === 验证 YZ ≈ I ===
I_est = Y * Z;               % 计算 YZ
I_ref = eye(size(Y));        % 理论单位矩阵
Delta = I_est - I_ref;          % 误差矩阵

% 误差指标
err_F = norm(Delta, 'fro');
err_max = max(abs(Delta(:)));
fprintf('[YZ-I] |·|_max = %.3e, ||·||_F = %.3e\n', err_max, err_F);

% === 导出为 CSV 文件 ===
outFile = fullfile(pwd, 'YZ_minus_I.csv');  % 保存路径
writematrix(full(Delta), outFile);

fprintf('已导出到: %s\n', outFile);
