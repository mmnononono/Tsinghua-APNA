function [V_eq, Z_eq] = Thevenin(system, port_list, V, Z)
define_constants;

bus = system.bus;
bus_no = bus(:, BUS_I);
nb = length(bus_no);
ng = nb + 1;

V = V(:);
Vpad = [V; 0];

Zpad = zeros(ng);
Zpad(1:nb,1:nb) = Z;

nport = size(port_list,1);
M = zeros(ng, nport);

for k = 1:nport
    a = port_list(k,1);
    b = port_list(k,2);
    if a ~= ng
        M(bus_no==a,k) = 1;
    else
        M(ng,k) = 1;
    end
    if b ~= ng
        M(bus_no==b,k) = -1;
    else
        M(ng,k) = -1;
    end
end

V_eq = M.' * Vpad;
Z_eq = M.' * Zpad * M;
end
