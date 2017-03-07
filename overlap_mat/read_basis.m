function [ mol_basis,nbasis_pg,nbasis,ncoef_pg] = read_basis( mol,natoms,basis_name )

nbasis_pg = 0;
nbasis = 0;
ncoef_pg =0;

for iatom = 1:natoms
  atom = mol(iatom).atom;
  
  fprintf('  %s basis functions \n',atom);
  if strcmp(basis_name,'sto3g')
    iatom_bs = sto3g(atom);
  elseif strcmp(basis_name,'ccpvdz')
    iatom_bs = ccpvdz(atom);
      elseif strcmp(basis_name,'ccpvdz_2')
    iatom_bs = ccpvdz_2(atom);
  elseif strcmp(basis_name,'augccpvdz')
    iatom_bs = augccpvdz(atom);
  else
    fprintf('\n %s is not implemented yet \n\n',basis_name);
    break
  end
  
  for ibs = 1:length(iatom_bs)
    n_ibs = numel(iatom_bs(ibs).expo);
    n_coef = size(iatom_bs(ibs).coef,2);
    
    if strcmp(iatom_bs(ibs).name,'P')
      % there are three p orbitals
      nbasis_pg = nbasis_pg+n_ibs*3;
      nbasis = nbasis+n_coef*3;
      ncoef_pg = ncoef_pg+n_ibs*3*n_coef;
    elseif strcmp(iatom_bs(ibs).name,'D')
      % there are six d orbitals
      nbasis_pg = nbasis_pg+n_ibs*6;
      nbasis = nbasis+n_coef*6;
      ncoef_pg = ncoef_pg+n_ibs*6*n_coef;
    else 
      nbasis_pg = nbasis_pg+n_ibs;
      nbasis = nbasis+n_coef;
      ncoef_pg = ncoef_pg+n_ibs*n_coef;
    end

    mol_basis(iatom).basis(ibs).name = iatom_bs(ibs).name;
    mol_basis(iatom).basis(ibs).expo = iatom_bs(ibs).expo;
    mol_basis(iatom).basis(ibs).coef = iatom_bs(ibs).coef;
  end
  
end


end

