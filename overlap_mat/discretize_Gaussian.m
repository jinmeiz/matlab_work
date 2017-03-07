function [ G_x,G_y,G_z ] = discretize_Gaussian( x,y,z,n,nbasis_pg,natoms,...
                                                mol_basis,Axyz )
% discretization of primitive Gaussian functions
% b: size of box
% nbasis_pg: # of primitive Gaussian functions
% natoms: # of atoms
% mol_basis: basis information
% Axyz: coordinates of the centers of functions on 3D

G_x = zeros(n,nbasis_pg);
G_y = zeros(n,nbasis_pg);
G_z = zeros(n,nbasis_pg);

idx_bs = 1;
for iatom = 1:natoms
  iatom_bs = mol_basis(iatom).basis;
  
  for ibs = 1:length(iatom_bs)
      
    ibs_name = iatom_bs(ibs).name;
    ibs_expo = iatom_bs(ibs).expo;
    ibs_nbs =length(ibs_expo);
    
    if strcmp(ibs_name,'P')
      % there are three p orbitals
      for i = 1:ibs_nbs
        for p_idx = 1:3
          G_x(:,idx_bs) = Gp_1D(x,Axyz(1,idx_bs),ibs_expo(i),p_idx-1);
          G_y(:,idx_bs) = Gp_1D(y,Axyz(2,idx_bs),ibs_expo(i),p_idx-2);
          G_z(:,idx_bs) = Gp_1D(z,Axyz(3,idx_bs),ibs_expo(i),p_idx-3);
          
          idx_bs = idx_bs+1;
        end
      end
    elseif strcmp(ibs_name,'D')
      % there are six d orbitals
      for i = 1:ibs_nbs
        % dxx 
        G_x(:,idx_bs) = Gd_1D_1(x,Axyz(1,idx_bs),ibs_expo(i),0);
        G_y(:,idx_bs) = Gd_1D_1(y,Axyz(2,idx_bs),ibs_expo(i),1);
        G_z(:,idx_bs) = Gd_1D_1(z,Axyz(3,idx_bs),ibs_expo(i),1);          
        idx_bs = idx_bs+1;
        % dxy
        G_x(:,idx_bs) = Gd_1D_2(x,Axyz(1,idx_bs),ibs_expo(i),0);
        G_y(:,idx_bs) = Gd_1D_2(y,Axyz(2,idx_bs),ibs_expo(i),0);
        G_z(:,idx_bs) = Gd_1D_2(z,Axyz(3,idx_bs),ibs_expo(i),1);
        idx_bs = idx_bs+1;
        % dxz
        G_x(:,idx_bs) = Gd_1D_2(x,Axyz(1,idx_bs),ibs_expo(i),0);
        G_y(:,idx_bs) = Gd_1D_2(y,Axyz(2,idx_bs),ibs_expo(i),1);
        G_z(:,idx_bs) = Gd_1D_2(z,Axyz(3,idx_bs),ibs_expo(i),0);
        idx_bs = idx_bs+1;
        
        % dyy
        G_x(:,idx_bs) = Gd_1D_1(x,Axyz(1,idx_bs),ibs_expo(i),1);
        G_y(:,idx_bs) = Gd_1D_1(y,Axyz(2,idx_bs),ibs_expo(i),0);
        G_z(:,idx_bs) = Gd_1D_1(z,Axyz(3,idx_bs),ibs_expo(i),1);        
        idx_bs = idx_bs+1;
        % dyz
        G_x(:,idx_bs) = Gd_1D_2(x,Axyz(1,idx_bs),ibs_expo(i),1);
        G_y(:,idx_bs) = Gd_1D_2(y,Axyz(2,idx_bs),ibs_expo(i),0);
        G_z(:,idx_bs) = Gd_1D_2(z,Axyz(3,idx_bs),ibs_expo(i),0);
        idx_bs = idx_bs+1;  
        
        % dzz
        G_x(:,idx_bs) = Gd_1D_1(x,Axyz(1,idx_bs),ibs_expo(i),1);
        G_y(:,idx_bs) = Gd_1D_1(y,Axyz(2,idx_bs),ibs_expo(i),1);
        G_z(:,idx_bs) = Gd_1D_1(z,Axyz(3,idx_bs),ibs_expo(i),0);       
        idx_bs = idx_bs+1; 
      end
    elseif strcmp(ibs_name,'S') 
      for i = 1:ibs_nbs
        G_x(:,idx_bs) = Gs_1D(x,Axyz(1,idx_bs),ibs_expo(i));
        G_y(:,idx_bs) = Gs_1D(y,Axyz(2,idx_bs),ibs_expo(i));
        G_z(:,idx_bs) = Gs_1D(z,Axyz(3,idx_bs),ibs_expo(i));
        
        idx_bs = idx_bs+1;
      end
    else
      disp([ibs_name ' is not implemented']);
    end
  
  end
end


end

