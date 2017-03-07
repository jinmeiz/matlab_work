function [ g_1d ] = Gp_1D(r_d,A_d,alpha,factor)
% 1D primitative p-type Gaussian function 
% (r_d-A_d) e^(-aplha*(r_d-A_d)^2)
% A_d: center, r_d: coordinates on the d dimention

N = (128*alpha^5/pi^3)^(1/4/3);

if factor == 0
  beta = 1;
else
  beta = 0;
end

g_1d = N ...
       * (r_d-A_d).^beta ...
      .* exp(-alpha*((r_d-A_d).^2));
end

