function Y_new = Compensation_Y(Y, f_node, t_node, r, x, b, method)

n = size(Y, 1);
Y_new = Y;

z  = r + 1j * x;
y  = 1 / z;
yb = 1j * b / 2;

switch method
    case 'branch'
        Y_new(f_node, f_node) = Y_new(f_node, f_node) + y + yb;
        Y_new(t_node, t_node) = Y_new(t_node, t_node) + y + yb;
        Y_new(f_node, t_node) = Y_new(f_node, t_node) - y;
        Y_new(t_node, f_node) = Y_new(t_node, f_node) - y;

    case 'node'
        e_f = zeros(n, 1); e_f(f_node) = 1;
        e_t = zeros(n, 1); e_t(t_node) = 1;

        A  = [e_f.'; e_t.'];
        Yb = [y + yb, -y; -y, y + yb];

        Y_new = Y + A.' * Yb * A;
end
end
