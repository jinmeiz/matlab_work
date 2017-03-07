function [ C_pg2cg ] = compute_coef_pg( natoms,nbasis,ncoef_pg,mol_basis )
% compute the coefficient matrix C for primitive Gaussian

C_pg2cg = zeros(ncoef_pg,nbasis);

s = 1;
idx_s = 1;
for iatom = 1:natoms
  iatom_bs = mol_basis(iatom).basis;

  for ibs = 1:length(iatom_bs)
    ibs_name = iatom_bs(ibs).name;
    ibs_coef_tol = iatom_bs(ibs).coef;
        
    n_coef = size(ibs_coef_tol,2);
    for icoef = 1:n_coef
      ibs_coef = ibs_coef_tol(:,icoef);
      ibs_nbs =length(ibs_coef);
      
      if strcmp(ibs_name,'S')
        C_pg2cg(idx_s:idx_s+ibs_nbs-1,s) = ibs_coef;
        s = s+1;
        idx_s = idx_s+ibs_nbs;

      elseif strcmp(ibs_name,'P') 
        ip_offset = 0;
        for ifunc = 1:ibs_nbs
          coef_ifunc = ibs_coef(ifunc);
          C_pg2cg(idx_s+ip_offset,s) = coef_ifunc;
          C_pg2cg(idx_s+1+ip_offset,s+1) = coef_ifunc;
          C_pg2cg(idx_s+2+ip_offset,s+2) = coef_ifunc;
          ip_offset = ifunc*3;
        end
        idx_s = idx_s+ibs_nbs*3;
        s = s+3;

      elseif strcmp(ibs_name,'D') 
        ip_offset = 0;
        for ifunc = 1:ibs_nbs
          coef_ifunc = ibs_coef(ifunc);
          C_pg2cg(idx_s+ip_offset,s) = coef_ifunc;
          C_pg2cg(idx_s+1+ip_offset,s+1) = coef_ifunc;
          C_pg2cg(idx_s+2+ip_offset,s+2) = coef_ifunc;
          C_pg2cg(idx_s+3+ip_offset,s+3) = coef_ifunc;
          C_pg2cg(idx_s+4+ip_offset,s+4) = coef_ifunc;
          C_pg2cg(idx_s+5+ip_offset,s+5) = coef_ifunc;
          ip_offset = ifunc*6;
        end
        idx_s = idx_s+ibs_nbs*6;
        s = s+6;
      else
      disp([ibs_name ' is not implemented']);
      end  % end if statement
      
    end % end of looping over num of coefficients

  end
end 


end

