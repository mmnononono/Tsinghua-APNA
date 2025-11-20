function Y_factor_new = Factorization_Modification(Y_factor, Y, f_node, t_node, r, x, b, method)

[L, D, U] = unpack_factor(Y_factor);

n = size(Y, 1);

switch method
    case 'ReFactorization'
        S1 = [f_node, t_node];
        S2 = Build_Path_Set(L, U, S1);
        I2 = min(S2):max(S2);

        Y22 = Y(I2, I2);
        fi = find(I2 == f_node, 1);
        ti = find(I2 == t_node, 1);

        Y22_new = Compensation_Y(Y22, fi, ti, r, x, b, 'branch');
        [L22, D22, U22] = My_makeLDU(Y22_new);

        L(I2,I2) = L22;
        D(I2,I2) = D22;
        U(I2,I2) = U22;

    case 'rank1'
        y  = 1/(r+1j*x);
        yb = 1j*b/2;

        LDU = {{f_node, t_node, y}, {f_node, f_node, yb}, {t_node, t_node, yb}};

        for k = 1:3
            [L, D, U] = LDURank1_Update(L, D, U, LDU{k}{1}, LDU{k}{2}, LDU{k}{3});
        end
end

Y_factor_new = pack_factor(L, D, U);
end


function S2 = Build_Path_Set(L, U, S1)
n = size(L, 1);
v = false(n,1);

q = unique(S1(:));
v(q) = true;
h=1;

while h <= numel(q)
    k = q(h); h=h+1;
    nei = unique([find(L(k,:)~=0), find(U(k,:)~=0), find(L(:,k)~=0).', find(U(:,k)~=0).']);
    for x = nei
        if ~v(x)
            v(x) = true;
            q(end+1,1) = x;
        end
    end
end
S2 = find(v);
end


function [L, D, U] = LDURank1_Update(L, D, U, i, j, a)

n = size(L,1);

M = zeros(n,1);
N = zeros(n,1);

if i~=j
    M(i)=1; M(j)=-1;
    N(i)=1; N(j)=-1;
else
    M(i)=1; N(i)=1;
end

for k = 1:n
    d = D(k,k);
    mk = M(k);
    nk = N(k);

    d2 = d + mk*a*nk;
    D(k,k) = d2;

    if k < n
        idx = k+1:n;

        l = L(idx,k);
        u = U(k,idx);

        N1 = N(idx) - nk*u.';
        U(k,idx) = u + (mk*a/d2)*N1.';

        M1 = M(idx) - l*mk;
        L(idx,k) = l + (a*nk/d2)*M1;

        a = a - a*nk*(1/d2)*mk*a;
        M(idx) = M1;
        N(idx) = N1;
    end

    M(k)=0; N(k)=0;
end
end
