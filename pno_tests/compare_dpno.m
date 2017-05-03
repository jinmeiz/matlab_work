%% compare the subspaces between different PNO coefficients
tic;

mol_name = 'h2o';
ndocc = 4;
nvir = 8;

%% construct PNO coefficients from files
fprintf('\n constructing PNO coefficients for %s\n', mol_name);

Dab_ij = zeros(nvir,nvir,ndocc,ndocc);
Dab_ij_new = zeros(nvir,nvir,ndocc,ndocc);
for iter = 1:3
  fprintf('\n%2dth iteration results: \n',iter-1);
  for i = 1:ndocc
  
    for j = 1:i
      f_name = strcat(strcat(strcat(strcat(strcat('C_es_', ...
                   int2str(i)),int2str(j)),'_'),int2str(iter-1)),'.out');
      f_new_name = strcat(strcat(strcat(strcat(strcat('C_es_', ...
                   int2str(i)),int2str(j)),'_'),int2str(iter-1)),'new.out');
      %fprintf('  reading %10s\n',f_name);
      Dab = load([strcat(strcat('./',mol_name),'/') f_name]); 
      n_pno = size(Dab,2);
      Dab_ij(:,1:n_pno,i,j) = Dab;
      %fprintf('  reading %10s\n',f_new_name);
      Dab_new = load([strcat(strcat('./',mol_name),'/') f_new_name]); 
      n_pno_new = size(Dab_new,2);
      Dab_ij_new(:,1:n_pno_new,i,j) = Dab_new;
      
      % compute the angle between two PNO coefficients
      theta = subspace(Dab,Dab_new);
      fprintf('  for ij (%d,%d) pair, rank: %d vs. %d, theta = %f\n',i,j,n_pno,n_pno_new,theta);
    end
  
  end

end

toc;

