clear; clc;
v = ver('MATLAB'); fprintf('MATLAB %s %s\n', v.Name, v.Version);
mpver;

%Y0
mpc = loadcase('case39')
[Y0, Y] = My_makeY(mpc);
[Ybus, ~, ~] = makeYbus(mpc);
D = Y - Ybus;
fprintf('[case39] |Y-Ybus|_max=%.3e, ||Y-Ybus||_F=%.3e\n', full(max(abs(D(:)))), norm(D,'fro'));

%tap=0.93
mpc2 = loadcase('case39'); mpc2.branch(6,9) = 0.93;
[Y0_mod, Y_mod] = My_makeY(mpc2);
[Ybus_mod, ~, ~] = makeYbus(mpc2);
Dmod = Y_mod - Ybus_mod;
fprintf('[case39|tap=0.93] |Y-Ybus|_max=%.3e, ||Y-Ybus||_F=%.3e\n', full(max(abs(Dmod(:)))), norm(Dmod,'fro'));

%3更大算例系统
cases   = {'case14','case30','case57','case118'};
nRepeat = 8;               
fprintf('%-8s %6s %12s %12s %10s %12s\n', 'case','nbus','|Δ|_max','||Δ||_F','t(My)s','t(Ybus)s');

for k = 1:numel(cases)
    cname = cases{k};
    mpc0  = loadcase(cname);    
    mpc   = ext2int(mpc0);     

    [~, Ymy] = My_makeY(mpc0);
    [Ybus,~,~] = makeYbus(mpc);

    D   = Ymy - Ybus;
    d1  = full(max(abs(D(:))));
    dF  = norm(full(D),'fro');

    t1 = 0; t2 = 0;
    for r = 1:nRepeat
        tic; [~,~] = My_makeY(mpc0); t1 = t1 + toc;
        tic; [~,~,~] = makeYbus(mpc); t2 = t2 + toc;
    end
    t1 = t1/nRepeat; t2 = t2/nRepeat;

    fprintf('%-8s %6d %12.3e %12.3e %10.6f %12.6f\n', cname, size(mpc.bus,1), d1, dF, t1, t2);
end