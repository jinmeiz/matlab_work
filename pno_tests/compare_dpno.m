%% compare the subspaces between different PNO coefficients
% there are six different PNO coefficients
tic;

%% calculation information
% H2O/6-31G calculaton: # of occ: 4, # of vir: 8
mol_name = 'h2o';
ndocc = 4;
nvir = 8;

%% construct PNO coefficients from files
fprintf('\n constructing PNO coefficients for %s\n', mol_name);

Dab_ij_tem1 = zeros(nvir,nvir,ndocc,ndocc);
Dab_ij_tem2 = zeros(nvir,nvir,ndocc,ndocc);

theta_1 = zeros(ndocc,ndocc);
theta_2 = zeros(ndocc,ndocc);
idx = 1;
for iter = 1:3
  fprintf('\n%2dth iteration results: \n',iter-1);
  
  Dab_ij_tem2(:,:) = 0;
  for i = 1:ndocc  
    for j = 1:i
      f_name = strcat(strcat(strcat(strcat(strcat('C_es_', ...
                   int2str(i)),int2str(j)),'_'),int2str(iter-1)),'.out');
      f_new_name = strcat(strcat(strcat(strcat(strcat('C_es_', ...
                   int2str(i)),int2str(j)),'_'),int2str(iter-1)),'new.out');
      %fprintf('  reading %10s\n',f_name);
      Dab = load([strcat(strcat('./',mol_name),'/') f_name]); 
      
      if (iter > 1)
        % compute the angle between previous and current PNO coefficients
        theta_2(i,j) = subspace(Dab_ij_tem1(:,:,i,j),Dab);
%         if (i == 4 && j == 1 && iter == 3)
%           disp(i);
%           disp(j);
%           disp(Dab);
%           disp(Dab_ij_tem1(:,:,i,j));
%           disp(theta_2(i,j))
%         end
      end
      
      %fprintf('  reading %10s\n',f_new_name);
      Dab_new = load([strcat(strcat('./',mol_name),'/') f_new_name]); 
      
      % compute the angle between two PNO coefficients
      theta_1(i,j) = subspace(Dab,Dab_new);
      
      % store the last PNO coefficients
      Dab_ij_tem2(:,1:size(Dab_new,2),i,j) = Dab_new;
      
%       if (i == 4 && j == 2 && iter == 3)
%         disp(f_name);
%         disp(f_new_name);
%         disp(i);
%         disp(j);
%         disp(idx);
%         disp(Dab);
%         disp(Dab_new);
%         %disp(Dab_ij(:,1:n_pno,i,j,iter));
%       end
    end  
  end
  Dab_ij_tem1 = Dab_ij_tem2;
  
  if (iter > 1)
  fprintf('\n%2dth comparison (cross iterations): \n',idx);
  for i = 1:ndocc  
    for j = 1:i
      fprintf('  for ij (%d,%d) pair, theta = %f\n',i,j,theta_2(i,j));
    end
  end
  idx = idx+1;
  end
  
  fprintf('\n%2dth comparison: \n',idx);
  for i = 1:ndocc  
    for j = 1:i
      fprintf('  for ij (%d,%d) pair, theta = %f\n',i,j,theta_1(i,j));
    end
  end
  idx = idx+1;

end

toc;

