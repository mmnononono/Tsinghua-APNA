function [V_eq_change, Z_eq_change] = Thevenin_Change(system, port_list, V, Z, f_node, t_node, V_eq, Z_eq)
define_constants;
bus    = system.bus;
branch = system.branch;
bus_no = bus(:, BUS_I);
nb     = length(bus_no);
ng     = nb + 1;

id = find((branch(:,F_BUS)==f_node & branch(:,T_BUS)==t_node) | (branch(:,F_BUS)==t_node & branch(:,T_BUS)==f_node), 1);

r = branch(id, BR_R);
x = branch(id, BR_X);
y = 1/(r + 1j*x);
Delta_y = -y;

m = zeros(nb+1,1);
m(bus_no==f_node)= 1;
m(bus_no==t_node)=-1;

Zpad = zeros(ng);
Zpad(1:nb,1:nb) = Z;

Z_alpha = Zpad * m;
Z_aa    = m.' * Z_alpha;
Y_aa    = 1 / (1/Delta_y + Z_aa);

nport = size(port_list,1);
M = zeros(ng,nport);
for k = 1:nport
    a = port_list(k,1);
    b = port_list(k,2);

    if a ~= ng, M(bus_no==a,k)= 1; else, M(ng,k)= 1; end
    if b ~= ng, M(bus_no==b,k)=-1; else, M(ng,k)=-1; end
end

Z_La = M.' * Z_alpha;

Vpad = [V(:); 0];
V_alpha_0 = m.' * Vpad;

V_eq_change = V_eq - Z_La * Y_aa * V_alpha_0;
Z_eq_change = Z_eq - Z_La * Y_aa * Z_La.';
end
