function Z = My_makeZ(system)
% My_makeZ — 基于"对地支路初始化 + 连支追加"的节点阻抗矩阵生成
% 说明：
%   - 若 system.EPS_G 存在则采用；否则默认 1e-9
%   - 若 system.coupled_groups 存在则做 rank-m Woodbury；否则逐条 rank-1

    % ===== 阻尼设置（可调） =====
   
        EPS_G = 1e-8;
   

    % ===== 索引常量 =====
    define_constants;  % BUS_I, GS, BS, F_BUS, T_BUS, BR_R, BR_X, BR_B, TAP, SHIFT, BR_STATUS

    baseMVA = system.baseMVA;
    bus     = system.bus;
    branch  = system.branch;

    nb = size(bus,1);
    nl = size(branch,1);

    % ===== (A) 初始化 Z^(0) = diag(1./y_sh) =====
    ysh = complex(zeros(nb,1));

    % 1) 母线并联导纳（与地）
    ysh = ysh + (bus(:,GS) + 1j*bus(:,BS)) / baseMVA;

    % 2) 线路两端充电电纳一半（各自对地）
    on  = branch(:,BR_STATUS) > 0;
    if any(on)
        f_on  = branch(on, F_BUS);
        t_on  = branch(on, T_BUS);
        Bb_on = branch(on, BR_B);
        y_half = 1j * Bb_on / 2;
        for k = 1:numel(y_half)
            ysh(f_on(k)) = ysh(f_on(k)) + y_half(k);
            ysh(t_on(k)) = ysh(t_on(k)) + y_half(k);
        end
    end

    % 防止零地导纳（数值稳健）
    ysh = ysh + EPS_G;

    Z = spdiags(1 ./ ysh, 0, nb, nb);  % 初始对角 Z

    % ===== (B) 追加"连支" =====
    coupled_groups = {};
    if isfield(system, 'coupled_groups') && ~isempty(system.coupled_groups)
        coupled_groups = system.coupled_groups;
    end

    % 标记耦合支路，避免重复
    coupled_mask = false(nl,1);
    for g = 1:numel(coupled_groups)
        members = coupled_groups{g}.members(:);
        assert(all(members>=1 & members<=nl), 'coupled members 越界');
        coupled_mask(members) = true;
    end

    % —— 无耦合：逐条 rank-1
    on_uncoupled = on & ~coupled_mask;
    if any(on_uncoupled)
        idx = find(on_uncoupled).';
        for k = idx
            z = complex(branch(k,BR_R), branch(k,BR_X));
            if ~isfinite(z) || z == 0, continue; end
            M = make_M_vector(nb, branch(k,:), F_BUS, T_BUS, TAP, SHIFT);

            ZM = Z * M;
            S  = z + (M.' * ZM);            % 标量
            if ~isfinite(S) || abs(S) == 0, continue; end
            Z  = Z - (ZM * (ZM.')) / S;
        end
    end

    % —— 耦合组：rank-m Woodbury
    for g = 1:numel(coupled_groups)
        grp = coupled_groups{g};
        members = grp.members(:).';
        Zb = grp.Zb;

        assert(all(on(members)), '耦合组包含未投运支路');
        m = numel(members);
        assert(all(size(Zb)==[m,m]), 'Zb 尺寸必须为 m×m');

        Mb = complex(zeros(nb,m));
        for j = 1:m
            Mb(:,j) = make_M_vector(nb, branch(members(j),:), F_BUS, T_BUS, TAP, SHIFT);
        end

        ZMb = Z * Mb;
        S   = Zb + (Mb.' * ZMb);
        X   = S \ (ZMb.');
        Z   = Z - ZMb * X;
    end
end

function M = make_M_vector(nb, br, F_BUS, T_BUS, TAP, SHIFT)
    f = br(F_BUS);  t = br(T_BUS);
    tap   = br(TAP);   if tap==0 || ~isfinite(tap), tap = 1; end
    shift = br(SHIFT); if ~isfinite(shift), shift = 0; end
    a = tap * exp(1j * deg2rad(shift));
    M    = complex(zeros(nb,1));
    M(f) =  1 / a;
    M(t) = -1;
end
