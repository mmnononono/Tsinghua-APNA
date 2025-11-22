clear; clc;
define_constants;

system = case39;
bus = system.bus;
bus_no = bus(:, BUS_I);
nbus = length(bus_no);

Bound = [1 3 26];
External = [2 25 30 37];
Internal = setdiff(bus_no.', [Bound External]);

Ib = [-0.1249+1j*0.08572; 0.0741-1j*0.4595; -0.0886+1j*0.1693];
Ie = [0.3969+1j*0.4041; 0; 0.2135+1j*0.4160; 0.1038-1j*0.8289];
Ii = zeros(length(Internal), 1);

[Yw, Iw, Ub_w, Ui_w] = Ward(system, External, Bound, Internal, Ie, Ib, Ii);
U_BI_w = Yw \ Iw;
err_W1 = norm([Ub_w; Ui_w] - U_BI_w, inf);

[~, Y] = My_makeY(system);
I_full = zeros(nbus,1);
I_full(ismember(bus_no, Bound)) = Ib;
I_full(ismember(bus_no, External)) = Ie;

U_full = Y \ I_full;

Ub_f = U_full(ismember(bus_no, Bound));
Ui_f = U_full(ismember(bus_no, Internal));

err_W2 = norm(Ub_w - Ub_f, inf);
err_W3 = norm(Ui_w - Ui_f, inf);

port_list = [13 14; 13 40; 14 40];

Z = My_makeZ(system);
[V_eq, Z_eq] = Thevenin(system, port_list, U_full, Z);

v13 = U_full(bus_no==13);
v14 = U_full(bus_no==14);
U_port_full = [v13-v14; v13; v14];

err_T1 = norm(V_eq - U_port_full, inf);

np = size(port_list,1);
ng = nbus + 1;
M_alpha = zeros(ng, np);
for k = 1:np
    a = port_list(k,1);
    b = port_list(k,2);
    if a ~= ng
        M_alpha(bus_no==a,k) = 1;
    else
        M_alpha(ng,k) = 1;
    end
    if b ~= ng
        M_alpha(bus_no==b,k) = -1;
    else
        M_alpha(ng,k) = -1;
    end
end

Zpad = zeros(ng);
Zpad(1:nbus,1:nbus) = Z;
err_T2 = norm(Z_eq - M_alpha.'*Zpad*M_alpha, inf);

[V_eq_change, Z_eq_change] = Thevenin_Change(system, port_list, U_full, Z, 4, 14, V_eq, Z_eq);

sys2 = system;
id = find((sys2.branch(:,F_BUS)==4 & sys2.branch(:,T_BUS)==14) | ...
          (sys2.branch(:,F_BUS)==14 & sys2.branch(:,T_BUS)==4), 1);
sys2.branch(id, BR_STATUS) = 0;

Z2 = My_makeZ(sys2);
I_inj = Z \ U_full;
U2 = Z2 * I_inj;

[V_eq2, Z_eq2] = Thevenin(sys2, port_list, U2, Z2);

err_C1 = norm(V_eq_change - V_eq2, inf);
err_C2 = norm(Z_eq_change - Z_eq2, inf);

fprintf("WARD 等值验证：\n");
fprintf("Ward 内部电压一致性：%.3e\n", err_W1);
fprintf("边界节点误差：%.3e\n", err_W2);
fprintf("内部节点误差：%.3e\n\n", err_W3);

fprintf("Thevenin 等值验证：\n");
fprintf("端口电压误差：%.3e\n", err_T1);
fprintf("端口阻抗误差：%.3e\n\n", err_T2);

fprintf("Thevenin 修正验证：\n");
fprintf("等值电压误差：%.3e\n", err_C1);
fprintf("等值阻抗误差：%.3e\n", err_C2);
