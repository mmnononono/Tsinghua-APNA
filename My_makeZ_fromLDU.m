function Z = My_makeZ_fromLDU(L, D, U)
  n = size(D, 1);
  Z = zeros(n, n);
  I = eye(n);
  d = diag(D);

  for j = 1:n
      e = I(:, j);
      y = zeros(n, 1);
      for k = 1:n
          s = e(k);
          for m = 1:k-1       
              s = s - L(k,m) * y(m);
          end
          y(k) = s;
      end
      w = y ./ d;              
      x = zeros(n, 1);
      for k = n:-1:1
          s = w(k);
          for m = k+1:n        
              s = s - U(k,m) * x(m);
          end
          x(k) = s;
      end
      Z(:, j) = x;            
  end
end
