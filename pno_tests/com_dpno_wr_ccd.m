%% compare the subspaces between different PNO coefficients
% there are six different PNO coefficients
% this script stores all the PNO coefficients and rank information
tic;

%% calculation information
tcut = '1e7';
%% H2O
mol_name = 'h2o';
ndocc = 4;
% % H2O/6-31G calculaton: # of occ: 4, # of vir: 8
% bs_name = '631g';
% nvir = 8;
% H2O/aug-cc-pVDZ calculaton: # of occ: 4, # of vir: 36
bs_name = 'augdz';
nvir = 36;
%% (H2O)2 
% mol_name = 'h2o_2';
% ndocc = 8;
% % (H2O)2/6-31G calculaton: # of occ: 8, # of vir: 16
% % bs_name = '631g';
% % nvir = 16;
% (H2O)2/aug-cc-pVDZ calculaton: # of occ: 8, # of vir: 72
% bs_name = 'augdz';
% nvir = 72;

%% construct PNO coefficients from files
fprintf('\n constructing PNO coefficients for %s %s with %s\n', ...
  mol_name,bs_name,tcut);

Dab_ij = zeros(nvir,nvir,ndocc,ndocc,6);
n_pno = zeros(ndocc,ndocc,6);

for iter = 1:6
 
  for i = 1:ndocc  
    for j = 1:i
      f_name = strcat(strcat(strcat(strcat(strcat('C_es_', ...
                   int2str(i)),int2str(j)),'_'),int2str(iter-1)),'th.out');
      %fprintf('  reading %10s\n',f_name);
      path_name = strcat(strcat(strcat(strcat(strcat(strcat('./',mol_name), ...
      '/'),'bs_'),bs_name),'_'),tcut);
      Dab = load([strcat(path_name,'/new/') f_name]); 
      n_pno(i,j,iter) = size(Dab,2);
           
      Dab_ij(:,1:n_pno(i,j,iter),i,j,iter) = Dab;  
    end  
  end

end

%% import CCD PNO coefficients
fprintf('\n constructing PNO coefficients from CCD amplitudes\n');
Dab_ij_ccd = zeros(nvir,nvir,ndocc,ndocc);
n_pno_ccd = zeros(ndocc,ndocc);

  for i = 1:ndocc  
    for j = 1:i
      f_name_ccd = strcat(strcat(strcat(strcat('C_es_', ...
                   int2str(i)),int2str(j)),'_'),'ccd.out');
      %fprintf('  reading %10s\n',f_name);
      path_name_ccd = strcat(strcat(strcat(strcat(strcat(strcat('./',mol_name), ...
      '/'),'bs_'),bs_name),'_'),tcut);
      Dab_ccd = load([strcat(path_name_ccd,'/ccd/') f_name_ccd]); 
      n_pno_ccd(i,j) = size(Dab_ccd,2);
           
      Dab_ij_ccd(:,1:n_pno_ccd(i,j),i,j) = Dab_ccd;  
    end  
  end
  
%% import MP2 PNO coefficients
fprintf('\n constructing PNO coefficients from MP2 amplitudes\n');
Dab_ij_mp2 = zeros(nvir,nvir,ndocc,ndocc);
n_pno_mp2 = zeros(ndocc,ndocc);

  for i = 1:ndocc  
    for j = 1:i
      f_name_mp2 = strcat(strcat(strcat(strcat('C_es_', ...
                   int2str(i)),int2str(j)),'_'),'mp2.out');
      %fprintf('  reading %10s\n',f_name);
      path_name_mp2 = strcat(strcat(strcat(strcat(strcat(strcat('./',mol_name), ...
      '/'),'bs_'),bs_name),'_'),tcut);
      Dab_mp2 = load([strcat(path_name_mp2,'/mp2/') f_name_mp2]); 
      n_pno_mp2(i,j) = size(Dab_mp2,2);
           
      Dab_ij_mp2(:,1:n_pno_mp2(i,j),i,j) = Dab_mp2;  
    end  
  end
  
%% compare MP2 and CCD PNOs
for i = 1:ndocc  
  for j = 1:i
    n_pno_1 = n_pno_mp2(i,j);
    n_pno_2 = n_pno_ccd(i,j);
    Dab_1 = reshape(Dab_ij_mp2(:,:,i,j),nvir,nvir);
    Dab_2 = reshape(Dab_ij_ccd(:,:,i,j),nvir,nvir);
    theta = subspace(Dab_1(:,1:n_pno_1),Dab_2(:,1:n_pno_2))/pi*180;
%     fprintf('  for ij (%d,%d) pair, rank: %2d vs. %2d, theta = %5.2f\n',i,j, ...
%       n_pno_2,n_pno_1,theta);
    fprintf('  %5.2f', theta);
  end
end
  
%% print out label
fprintf('\n%s','label');
for i = 1:ndocc  
  for j = 1:i
    fprintf('  (%d,%d)',i,j);
  end
end
      
%% compute the subspace between different PNO coefficients
idx = 1;
for iter = 1:6
 fprintf('\n%2d ',idx);
 
  for i = 1:ndocc  
    for j = 1:i
      n_pno_1 = n_pno(i,j,iter);
      n_pno_2 = n_pno_ccd(i,j);
      Dab_1 = reshape(Dab_ij(:,:,i,j,iter),nvir,nvir);
      Dab_2 = reshape(Dab_ij_ccd(:,:,i,j),nvir,nvir);
      theta = subspace(Dab_1(:,1:n_pno_1),Dab_2(:,1:n_pno_2))/pi*180;
%       fprintf('  for ij (%d,%d) pair, rank: %2d vs. %2d, theta = %5.2f\n',i,j, ...
%         n_pno_2,n_pno_1,theta);
       fprintf('  %5.2f', theta);
    end
  end
  
  idx = idx+1;
end

fprintf('\n');
toc;
