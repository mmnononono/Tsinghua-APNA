function [Y0, Y] = My_makeY(system)
    define_constants;
    system = ext2int(system);
    baseMVA = system.baseMVA;
    bus     = system.bus;
    branch  = system.branch;

    n = size(bus,1);
    on = find(branch(:, BR_STATUS) ~= 0);  
    N  = n + 1;                           
    g  = N;                               
    Y0 = sparse(N, N);

    for k = 1:length(on)
        f = branch(on(k), F_BUS);
        t = branch(on(k), T_BUS);
        r = branch(on(k), BR_R);
        x = branch(on(k), BR_X);
        bch = branch(on(k), BR_B);
        tap = branch(on(k), TAP);
        shift = branch(on(k), SHIFT);  

        if tap == 0
            tap = 1;
        end

        a = tap * exp(1j * deg2rad(shift));
        y = 1 / (r + 1j*x);

        Y0(f,f) = Y0(f,f) + y / abs(a)^2;
        Y0(t,t) = Y0(t,t) + y;
        Y0(f,t) = Y0(f,t) - y / conj(a);
        Y0(t,f) = Y0(t,f) - y / a;

        % 充电电纳对地
        y_end_f = 1j * (bch/2) / abs(a)^2;   
        y_end_t = 1j * (bch/2);              

        Y0(f,f) = Y0(f,f) + y_end_f;
        Y0(g,g) = Y0(g,g) + y_end_f;
        Y0(f,g) = Y0(f,g) - y_end_f;
        Y0(g,f) = Y0(g,f) - y_end_f;

        
        Y0(t,t) = Y0(t,t) + y_end_t;
        Y0(g,g) = Y0(g,g) + y_end_t;
        Y0(t,g) = Y0(t,g) - y_end_t;
        Y0(g,t) = Y0(g,t) - y_end_t;
    end

    % 母线GS, BS
    for i = 1:n
        y_sh = (bus(i, GS) + 1j * bus(i, BS)) / baseMVA;
        if y_sh ~= 0
            Y0(i,i) = Y0(i,i) + y_sh;
            Y0(g,g) = Y0(g,g) + y_sh;
            Y0(i,g) = Y0(i,g) - y_sh;
            Y0(g,i) = Y0(g,i) - y_sh;
        end
    end

    % Y
    Y = Y0(1:n, 1:n);
end
