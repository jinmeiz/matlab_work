function [ Spg_tol ] = compute_Spg_tol( natoms,mol_basis,Spg )
% Copy the missing parts of overlap matrix
% There are primitive Gaussians that have same exponents but different
% coefficients. Only the unique primitive functions are discretized.To
% transfer from primitive overlap to contracted overlap, the missing parts
% have to be put back.

idx_npg = 1;
idx_ins = 1;
for iatom = 1:natoms
  iatom_bs = mol_basis(iatom).basis;

  for ibs = 1:length(iatom_bs)
    ibs_name = iatom_bs(ibs).name;
    ibs_nbs = size(iatom_bs(ibs).coef,1);       
    ibs_ncoef = size(iatom_bs(ibs).coef,2);
    
    if strcmp(ibs_name,'S')
      idx_len = ibs_nbs;
    elseif strcmp(ibs_name,'P') 
      idx_len = ibs_nbs*3;
    elseif strcmp(ibs_name,'D')
      idx_len = ibs_nbs*6;
    else
      err_meg = strcat(ibs_name,' is not implemented');
      error(err_meg);
    end
    
    if abs(ibs_ncoef-1)>1e-06
      ins_entry(idx_ins).idx = idx_npg;
      ins_entry(idx_ins).len = idx_len;
      ins_entry(idx_ins).col = Spg(:,idx_npg:idx_npg+idx_len-1)*(ibs_ncoef-1);
      idx_ins = idx_ins+1;
    end
    
    idx_npg = idx_npg+idx_len;  
  end
end 

n_ins = length(ins_entry);

% add columns
inc_ins = 0;
Spg_tol = Spg;
for i = 1:n_ins
  idx_ins = ins_entry(i).idx+ins_entry(i).len-1+inc_ins;
  % add a column
  Spg_tol = [Spg_tol(:,1:idx_ins),ins_entry(i).col, ...
             Spg_tol(:,idx_ins+1:end)];
  
  inc_ins = inc_ins+ins_entry(i).len;
end

% add rows
inc_ins = 0;
for i = 1:n_ins
  idx_b = ins_entry(i).idx+inc_ins;
  idx_e = idx_b+ins_entry(i).len-1;
  % add a row
  ins_row = Spg_tol(idx_b:idx_e,:);
  Spg_tol = [Spg_tol(1:idx_e,:);ins_row;Spg_tol(idx_e+1:end,:)];
 
  inc_ins = inc_ins+ins_entry(i).len;
end

end

