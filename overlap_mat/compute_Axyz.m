function [ Axyz ] = compute_Axyz( mol,mol_basis,natoms,nbasis_pg )
% compute coordinates for Gaussian functions

Axyz = zeros(3,nbasis_pg);
for i = 1:3
  
  A_idx = 1;
  for iatom = 1:natoms
    iatom_bs = mol_basis(iatom).basis;
    
    nbs = 0;
    for ibs = 1:length(iatom_bs)
      n_ibs = numel(iatom_bs(ibs).expo);
      
      if strcmp(iatom_bs(ibs).name,'P')
        % there are three p orbitals
        nbs = nbs+n_ibs*3;
      elseif strcmp(iatom_bs(ibs).name,'D')
        % there are six d orbitals
        nbs = nbs+n_ibs*6;  
      else 
        nbs = nbs+n_ibs;
      end
    end
    
    Axyz(i,A_idx:A_idx+nbs-1) = mol(iatom).xyz(i);
    A_idx = A_idx+nbs;
  end
end


end

