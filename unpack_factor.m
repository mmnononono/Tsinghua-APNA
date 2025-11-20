function [L, D, U] = unpack_factor(F)
n = size(F,1);
L = eye(n);
U = eye(n);
D = diag(diag(F));
L(tril(true(n),-1)) = F(tril(true(n),-1));
U(triu(true(n),1))  = F(triu(true(n),1));
end
