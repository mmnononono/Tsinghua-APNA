function F = pack_factor(L, D, U)
n = size(D,1);
F = zeros(n);
F(1:n+1:end) = diag(D);
F(triu(true(n),1)) = U(triu(true(n),1));
F(tril(true(n),-1)) = L(tril(true(n),-1));
end

