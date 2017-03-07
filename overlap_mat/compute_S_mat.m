% compute the overlap matrix of contracted Gaussian functions
% use h2o as example

%% begining of input
% import molecule geometry and basis
mol_name = 'h2o';
basis_name = 'sto3g';
ndocc = 5;
fprintf('\n reading %s geometry \n',mol_name);
eval(strcat(mol_name,'_xyz')); 
%%% end of input

%% read basis functions
natoms = length(mol);
fprintf('\n reading basis functions: %s \n',basis_name);
[ mol_basis,nbasis_pg,nbasis,ncoef_pg ] = read_basis( mol,natoms,basis_name );

fprintf(' number of primitive functions: %d\n',nbasis_pg);
fprintf(' number of contracted functions: %d\n',nbasis);
fprintf(' number of coefficients for primitive functions: %d\n',ncoef_pg);

%% compute coordinates for Gaussian functions
Axyz = compute_Axyz( mol,mol_basis,natoms,nbasis_pg );

%% compute discretizated basis functions
% box size [-b, b]
b = 12.5;
% grid number and size
p = 12;
n = 2^p;
h = 2*b/n;
x = (-b:h:b)';
y = x;
z = x; 
fprintf('\n box size: [-%d,%d]\n grid number: %d \n',b,b,n);

% discretization of Gaussian functions
fprintf('\n discretizating basis functions \n');
[G_x,G_y,G_z] = discretize_Gaussian( x,y,z,n+1,nbasis_pg,natoms,mol_basis,Axyz );

% basis function products 
GG_x = zeros(n+1,nbasis_pg,nbasis_pg);
GG_y = zeros(n+1,nbasis_pg,nbasis_pg);
GG_z = zeros(n+1,nbasis_pg,nbasis_pg);
for i = 1:1:nbasis_pg
  for j = 1:1:nbasis_pg
    GG_x(:,i,j) = G_x(:,i).*G_x(:,j);
    GG_y(:,i,j) = G_y(:,i).*G_y(:,j);
    GG_z(:,i,j) = G_z(:,i).*G_z(:,j);
  end
end

tic;

%% compute coefficients for tansforming primitive Gaussians to contracted
% Gaussian
C_pg2cg = compute_coef_pg( natoms,nbasis,ncoef_pg,mol_basis );  

%% test: calculate the overlap matrix
fprintf('\n computing overlap matrix: \n');
fprintf('  Spq with primitive functions \n');
Spq = zeros(nbasis_pg,nbasis_pg);
for i = 1:1:nbasis_pg
  for j = 1:1:nbasis_pg
    Spq(i,j) =  trapz(x,GG_x(:,i,j)) ...
              * trapz(y,GG_y(:,i,j)) ...
              * trapz(z,GG_z(:,i,j));
  end
end

% test: compute overlap matrix with contracted gaussians
fprintf('  Spq with contracted functions \n');
if abs(ncoef_pg-nbasis_pg)>1e-6
  Spg_tol = compute_Spg_tol( natoms,mol_basis,Spq );
else
  Spg_tol = Spq;
end
Spq_cg = C_pg2cg'*Spg_tol*C_pg2cg

toc;