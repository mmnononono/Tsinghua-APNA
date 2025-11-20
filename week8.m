clear; clc;

define_constants;
mpc = case39;

[~, Y] = My_makeY(mpc);

r1 = 0.052;  x1 = 0.025;  b1 = 0.258;
r2 = 0.061;  x2 = 0.028;  b2 = 0.205;

Y1_branch = Compensation_Y(Y, 3, 15, r1, x1, b1, 'branch');
Y1_node   = Compensation_Y(Y, 3, 15, r1, x1, b1, 'node');

Y2_branch = Compensation_Y(Y1_branch, 2, 28, r2, x2, b2, 'branch');
Y2_node   = Compensation_Y(Y1_node,   2, 28, r2, x2, b2, 'node');

fprintf('Y1: ||branch - node||_inf = %.3e\n',  norm(Y1_branch - Y1_node, 'inf'));
fprintf('Y2: ||branch - node||_inf = %.3e\n\n', norm(Y2_branch - Y2_node, 'inf'));

[L0, D0, U0] = My_makeLDU(Y);
Y_factor0 = pack_factor(L0, D0, U0);

Y_factor_r1_rank1 = Factorization_Modification(Y_factor0, Y, 3, 15, r1, x1, b1, 'rank1');
Y_factor_r1_refac = Factorization_Modification(Y_factor0, Y, 3, 15, r1, x1, b1, 'ReFactorization');

Y_factor_r2_rank1 = Factorization_Modification(Y_factor_r1_rank1, Y1_branch, 2, 28, r2, x2, b2, 'rank1');
Y_factor_r2_refac = Factorization_Modification(Y_factor_r1_refac, Y1_branch, 2, 28, r2, x2, b2, 'ReFactorization');

[L1_r1, D1_r1, U1_r1] = unpack_factor(Y_factor_r1_rank1);
[L2_r1, D2_r1, U2_r1] = unpack_factor(Y_factor_r2_rank1);

[L1_rf, D1_rf, U1_rf] = unpack_factor(Y_factor_r1_refac);
[L2_rf, D2_rf, U2_rf] = unpack_factor(Y_factor_r2_refac);

Y1 = Y1_branch;
Y2 = Y2_branch;

fprintf('rank1   : Y1  ||LDU - Y||_inf = %.3e\n', norm(L1_r1*D1_r1*U1_r1 - Y1, 'inf'));
fprintf('rank1   : Y2  ||LDU - Y||_inf = %.3e\n', norm(L2_r1*D2_r1*U2_r1 - Y2, 'inf'));
fprintf('Refactor: Y1  ||LDU - Y||_inf = %.3e\n', norm(L1_rf*D1_rf*U1_rf - Y1, 'inf'));
fprintf('Refactor: Y2  ||LDU - Y||_inf = %.3e\n\n', norm(L2_rf*D2_rf*U2_rf - Y2, 'inf'));

Z  = My_makeZ_fromLDU(L0, D0, U0);
[L1, D1, U1] = My_makeLDU(Y1);
Zp = My_makeZ_fromLDU(L1, D1, U1);
[L2, D2, U2] = My_makeLDU(Y2);
Zpp = My_makeZ_fromLDU(L2, D2, U2);

fprintf('原始  ||YZ-I||_inf = %.3e\n',  norm(Y  * Z   - eye(size(Y)), 'inf'));
fprintf(' Y1,Z1 ||YZ-I||_inf = %.3e\n',  norm(Y1 * Zp  - eye(size(Y)), 'inf'));
fprintf('Y2,Z2 ||YZ-I||_inf = %.3e\n\n', norm(Y2 * Zpp - eye(size(Y)), 'inf'));

n = size(Y,1);

z1  = r1 + 1j*x1; y1  = 1 / z1; yb1 = 1j * b1 / 2;
e3  = zeros(n,1); e3(3)  = 1;
e15 = zeros(n,1); e15(15) = 1;

m1_1 = e3 - e15;  a1_1 = y1;
m1_2 = e3;        a1_2 = yb1;
m1_3 = e15;       a1_3 = yb1;

Z_aux_p = rank1_update(Z, m1_1, a1_1);
Z_aux_p = rank1_update(Z_aux_p, m1_2, a1_2);
Z_aux_p = rank1_update(Z_aux_p, m1_3, a1_3);

z2  = r2 + 1j*x2; y2  = 1 / z2; yb2 = 1j * b2 / 2;
e2  = zeros(n,1); e2(2)  = 1;
e28 = zeros(n,1); e28(28) = 1;

m2_1 = e2 - e28; a2_1 = y2;
m2_2 = e2;       a2_2 = yb2;
m2_3 = e28;      a2_3 = yb2;

Z_aux_pp = rank1_update(Z_aux_p,  m2_1, a2_1);
Z_aux_pp = rank1_update(Z_aux_pp, m2_2, a2_2);
Z_aux_pp = rank1_update(Z_aux_pp, m2_3, a2_3);

fprintf('矩阵求逆辅助定理验证:\n');
fprintf(' ||Z1 - Z_aux1||_inf = %.3e\n', norm(Zp  - Z_aux_p,  'inf'));
fprintf(' ||Z2 - Z_aux2||_inf = %.3e\n', norm(Zpp - Z_aux_pp, 'inf'));

function Z_new = rank1_update(Z, M, a)
    x = Z * M;
    alpha = 1/a + M' * x;
    Z_new = Z - x * ((1/alpha) * x.');
end
