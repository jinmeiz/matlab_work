function [ L ] = compute_trunc_CD_mat( A,err_target)
% truncated Cholesky factorization of matrix A
% err_target is the threshold
n = size(A,1);

% compute diagonal elements
d = diag(A);

m = 1;
err = abs(d(1));
npi = 1:n;
L = zeros(n,n);
while err > err_target
  % find i = aregmax(d(pi_j):j=m,m+1,...,N)
  [~,i] = max(d(npi(m:n)));
  temp = npi(m);
  npi(m) = npi(i+m-1);
  npi(i+m-1) = temp;
  
  % diagonal elements: Lrr
  L(m,npi(m)) = sqrt(d(npi(m)));
  
  for i = m+1:n
    L(m,npi(i)) = A(npi(i),npi(m))/L(m,npi(m));
    if m > 1
      L(m,npi(i)) = L(m,npi(i)) ...
                   -dot(L(1:m-1,npi(m)),L(1:m-1,npi(i)))/L(m,npi(m));
    end
    
    d(npi(i)) = d(npi(i))-L(m,npi(i))^2;
  end
  err = sum(d(npi(m+1:n)));
  
  m = m+1;
end

% truncate L dimension to m * N
L(~any(L,2),:) = [];

end