function [L,D,U]=My_makeLDU(Y)
n=size(Y,1);
A=Y;
for p=1:n-1
    for j=p+1:n
        A(p,j)=A(p,j)/A(p,p);
    end
    for i=p+1:n
        for j=p+1:n
        A(i,j)=A(i,j)-A(i,p)*A(p,j);
        end
    end
end
L=tril(A);
U=triu(A,1)+eye(n);
D = diag(diag(L));         
L = L / D;                 