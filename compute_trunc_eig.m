function [ v,d] = compute_trunc_eig( mat,t_value,rel )
% compute truncated eigenvalues (v) and eigenvectors (d) for a matrix
% eigenvalues below t_value will be truncated
% if 'rel' is true, relative theshold is used
% if 'rel' is true, absolute theshold is used

[v,d] = eig(mat);

d_diag = abs(diag(d));
if ~rel
  idx = d_diag>t_value;
else 
  max_ddiag = max(d_diag);
  idx = d_diag>max_ddiag*t_value;
end

v = v(:,idx);
d = d(idx,idx);

%% test
% mat_test = v*d/v;
% norm_err = norm(mat_test-mat);
% fprintf('norm error: %3.2E   truncated rank: %3d \n', norm_err,size(v,2));

end

