function [ p_1d ] = compute_p_1d( b,n,C0,M )
% compute the tensor approximation of 1/r on one dimension using sinc-quadrature
% box size: [-b, b]
% n: number of grids
% C0 & M: sinc-quadrature parameter (e.g. C0=1.6, M=55)

eta = C0*log(M)/M;
rank = M+1;

%% use improved sinc-quadrature
krange1 = 0:1:M;
krange2 = 1:1:M;
g_k = zeros(1,M+1);
g_k(1) = eta;
g_k(2:end) = 2*eta.*cosh(krange2.*eta);
ak = pi^(-1/2).*g_k;
t_k = sinh(krange1.*eta);

p_1d = zeros(n+1,rank);
% Nystrom scheme: point value is used
h = 2*b/n;
x = (-b:h:b)';
for m = 1:rank  
  % use point value
  p_1d(:,m) = ak(m).^(1/3).*exp(-x.^2.*t_k(m).^2)*h;
end
s
%% use collocation scheme (piece-wise function)
% p_1d = zeros(n,rank);
% i = 1:n+1;
% h = 2*b/n;
% x_a = -b+h.*(i-1);
% x_a = x_a(:);
% erfmat = erf(x_a*t_k(2:end));
% p_1d(:,1) = ak(1)^(1/3)*h;
% p_1d(:,2:M+1) = repmat(ak(2:end),n,1).^(1/3) ...
%                .*sqrt(pi)/2./repmat(t_k(2:end),n,1) ...
%                .*(erfmat(2:end,:)-erfmat(1:end-1,:));
end

