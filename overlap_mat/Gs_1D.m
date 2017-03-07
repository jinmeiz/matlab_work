% 1D primitative s-type Gaussian function 
% (r_d-A_d)^beta e^(-aplha*(r_d-A_d)^2)
function [ g_1d ] = Gs_1D(r_d,A_d,alpha)

%A_d: center, r_d: coordinates on the d dimention
N = (2*alpha/pi)^(1/4);
%N = 1;
g_1d = N * exp(-alpha*((r_d-A_d).^2));

end


