%% compare the subspaces between different PNO coefficients
% there are six different PNO coefficients
% this script stores all the PNO coefficients and rank information
tic;

%% calculation information
mol_name = 'h2o';
ndocc = 4;
tcut = '1e7';
% H2O/6-31G calculaton: # of occ: 4, # of vir: 8
% bs_name = '631g';
% nvir = 8;
% H2O/aug-cc-pVDZ calculaton: # of occ: 4, # of vir: 36
bs_name = 'augdz';
nvir = 36;

%% construct PNO coefficients from files
fprintf('\n constructing PNO coefficients for %s %s with %s\n', ...
  mol_name,bs_name,tcut);

Dab_ij = zeros(nvir,nvir,ndocc,ndocc,3);
Dab_ij_new = zeros(nvir,nvir,ndocc,ndocc,3);
n_pno = zeros(ndocc,ndocc,3);
n_pno_new = zeros(ndocc,ndocc,3);

for iter = 1:3
 
  for i = 1:ndocc  
    for j = 1:i
      f_name = strcat(strcat(strcat(strcat(strcat('C_es_', ...
                   int2str(i)),int2str(j)),'_'),int2str(iter-1)),'.out');
      f_new_name = strcat(strcat(strcat(strcat(strcat('C_es_', ...
                   int2str(i)),int2str(j)),'_'),int2str(iter-1)),'new.out');
      %fprintf('  reading %10s\n',f_name);
      path_name = strcat(strcat(strcat(strcat(strcat(strcat('./',mol_name),'/'),'bs_'),bs_name),'_'),tcut);
      Dab = load([strcat(path_name,'/') f_name]); 
      n_pno(i,j,iter) = size(Dab,2);
      
      %fprintf('  reading %10s\n',f_new_name);
      Dab_new = load([strcat(path_name,'/') f_new_name]); 
      n_pno_new(i,j,iter) = size(Dab_new,2);
           
      Dab_ij(:,1:n_pno(i,j,iter),i,j,iter) = Dab;  
      Dab_ij_new(:,1:n_pno_new(i,j,iter),i,j,iter) = Dab_new;
    end  
  end

end

%% compute the subspace between different PNO coefficients
idx = 1;
for iter = 1:3
 
if (iter > 1)
  fprintf('\n%2dth comparison: \n',idx);
  
  for i = 1:ndocc  
    for j = 1:i
      n_pno_1 = n_pno(i,j,iter);
      n_pno_2 = n_pno_new(i,j,iter-1);
      Dab_1 = reshape(Dab_ij(:,:,i,j,iter),nvir,nvir);
      Dab_2 = reshape(Dab_ij_new(:,:,i,j,iter-1),nvir,nvir);
      fprintf('  for ij (%d,%d) pair, rank: %d vs. %d, theta = %e\n',i,j, ...
        n_pno_2,n_pno_1,subspace(Dab_1(:,1:n_pno_1),Dab_2(:,1:n_pno_2)));
    end
  end
  idx = idx+1;
end

  fprintf('\n%2dth comparison: \n',idx);
  for i = 1:ndocc  
    for j = 1:i
      n_pno_1 = n_pno(i,j,iter);
      n_pno_2 = n_pno_new(i,j,iter);
      Dab_1 = reshape(Dab_ij(:,:,i,j,iter),nvir,nvir);
      Dab_2 = reshape(Dab_ij_new(:,:,i,j,iter),nvir,nvir);
      fprintf('  for ij (%d,%d) pair, rank: %d vs. %d, theta = %e\n',i,j, ...
        n_pno_1,n_pno_2,subspace(Dab_1(:,1:n_pno_1),Dab_2(:,1:n_pno_2)));
    end
  end
  idx = idx+1;
end


% Dab_1 = Dab_ij(:,1:3,4,2,3);
% Dab_2 = Dab_ij_new(:,1:3,4,2,2);
% disp(Dab_1);
% disp(Dab_2);
% fprintf('  for ij (%d,%d) pair, theta = %f\n',4,2,subspace(Dab_2,Dab_1));

%% test
% Dab_1 = Dab_ij(:,:,4,1,3);
% Dab_2 = Dab_ij_new(:,:,4,1,2);
% disp(Dab_1);
% disp(Dab_2);
% fprintf('  for ij (%d,%d) pair, theta = %f\n',4,1,subspace(Dab_2,Dab_1));

toc;
