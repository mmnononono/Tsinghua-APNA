function Z = My_makeZ(system)

BUS_I=1; GS=5; BS=6; F_BUS=1; T_BUS=2; BR_R=3; BR_X=4; BR_B=5; TAP=9; SHIFT=10; BR_STATUS=11;

baseMVA = system.baseMVA;
bus     = system.bus;
branch  = system.branch;
nb = size(bus,1); nl = size(branch,1); ncol = size(branch,2);

tap = ones(nl,1);
if ncol >= TAP
    t = branch(:,TAP); nz = isfinite(t) & (t~=0); tap(nz)=t(nz);
end
shift_deg = zeros(nl,1);
if ncol >= SHIFT
    s = branch(:,SHIFT); s(isnan(s))=0; shift_deg=s;
end
on = true(nl,1);
if ncol >= BR_STATUS
    on = branch(:,BR_STATUS)>0;
end

%接地支路
ground_bus_id = max(bus(:,BUS_I)) + 1;
added_rows = {};

% Gs/Bs
if size(bus,2) >= max(GS,BS)
    for i=1:nb
        Yg=(bus(i,GS)+1i*bus(i,BS))/baseMVA;
        if abs(Yg)>0
            Zg=1/Yg; r=zeros(1,ncol);
            r([F_BUS T_BUS BR_R BR_X])=[bus(i,BUS_I) ground_bus_id real(Zg) imag(Zg)];
            if ncol>=BR_STATUS, r(BR_STATUS)=1; end
            added_rows{end+1}=r;
        end
    end
end

% 充电电纳 B/2
if size(branch,2) >= BR_B
    for k=1:nl
        if branch(k,BR_B)~=0
            Yc=1i*branch(k,BR_B)/2; Zc=1/Yc;
            rf=zeros(1,ncol); rf([F_BUS T_BUS BR_R BR_X])=[branch(k,F_BUS) ground_bus_id real(Zc) imag(Zc)];
            if ncol>=BR_STATUS, rf(BR_STATUS)=1; end
            added_rows{end+1}=rf;
            rt=rf; rt(F_BUS)=branch(k,T_BUS); added_rows{end+1}=rt;
        end
    end
end

branchE = branch;
if ~isempty(added_rows)
    branchE = [vertcat(added_rows{:}); branch];
end

f=int32(branchE(:,F_BUS)); t=int32(branchE(:,T_BUS));
zser=complex(branchE(:,BR_R),branchE(:,BR_X));
tapE=[ones(size(branchE,1)-nl,1); tap];
shiftE=[zeros(size(branchE,1)-nl,1); shift_deg];
aE=tapE.*exp(1i*pi/180*shiftE);
onE=[true(size(branchE,1)-nl,1); on];

%初始化 
Z=complex(0);
node_pos=containers.Map('KeyType','int32','ValueType','int32');
node_pos(ground_bus_id)=1;
applied=false(size(branchE,1),1);

    function add_tree_edge(k, known, newb, a)
        n=size(Z,1); ZN=complex(zeros(n+1,n+1)); ZN(1:n,1:n)=Z;
        ki=node_pos(known);
        if known==f(k)
            ZN(1:n,n+1)=Z(1:n,ki)/a; ZN(n+1,1:n)=Z(ki,1:n)/a;
            ZN(n+1,n+1)=zser(k)+Z(ki,ki)/(a^2);
        else
            ZN(1:n,n+1)=Z(1:n,ki)*a; ZN(n+1,1:n)=Z(ki,1:n)*a;
            ZN(n+1,n+1)=(a^2)*(Z(ki,ki)+zser(k));
        end
        Z=ZN; node_pos(newb)=n+1;
    end

    function add_link_edge(k,a)
        fi=node_pos(f(k)); tj=node_pos(t(k));
        M=complex(zeros(size(Z,1),1)); M(fi)=1/a; M(tj)=-1;
        v=Z*M; S=zser(k)+(M.'*v); Z=Z-(v*v.')/S;
    end

% 树支 
changed=true;
while changed
    changed=false;
    for k=1:size(branchE,1)
        if applied(k)||~onE(k), continue; end
        i=f(k); j=t(k);
        if xor(i==ground_bus_id,j==ground_bus_id)
            other=i; if i==ground_bus_id, other=j; end
            if ~isKey(node_pos,other)
                add_tree_edge(k,ground_bus_id,other,aE(k));
                applied(k)=true; changed=true;
            end
        end
    end
end

%  扩展树
changed=true;
while changed
    changed=false;
    for k=1:size(branchE,1)
        if applied(k)||~onE(k), continue; end
        i=f(k); j=t(k); ci=isKey(node_pos,i); cj=isKey(node_pos,j);
        if xor(ci,cj)
            known=i; newb=j; if ~ci&&cj, known=j; newb=i; end
            add_tree_edge(k,known,newb,aE(k));
            applied(k)=true; changed=true;
        end
    end
end

% 连支 
for k=1:size(branchE,1)
    if applied(k)||~onE(k), continue; end
    if isKey(node_pos,f(k))&&isKey(node_pos,t(k))
        add_link_edge(k,aE(k)); applied(k)=true;
    end
end

% 重排输出 
ids=int32(system.bus(:,BUS_I)); n0=numel(ids);
Z_out=complex(zeros(n0,n0));
for p=1:n0
    ip=node_pos(ids(p));
    for q=1:n0
        iq=node_pos(ids(q)); Z_out(p,q)=Z(ip,iq);
    end
end
Z=Z_out;
end
