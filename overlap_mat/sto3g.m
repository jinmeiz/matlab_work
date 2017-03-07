function [ basis ] = sto3g( atom )
% STO-3G basis for atom: H, C, N, and O
% AO basis functions are organized into shells: S, P, D, ...
% Each shell contains: 2l+1 spherical harmonic functions 
%                      (l = 0, 1, 2 for S, P, D) 
% For D shell, there are six cartesian AOs.

switch atom
  case 'H'
    % name of the shell
    basis(1).name = 'S';
    % exponents for AOs in the shell
    basis(1).expo = [3.42525091
                     0.62391373
                     0.16885540];
    % contraction coefficients for AOs in the shell
    basis(1).coef = [0.15432897
                     0.53532814
                     0.44463454];
  case 'C'
    basis(1).name = 'S';
    basis(1).expo = [71.616837000
                     13.045096000
                     3.530512200];
    basis(1).coef = [0.15432897
                     0.53532814
                     0.44463454];
    basis(2).name = 'S';
    basis(2).expo = [2.941249400
                     0.683483100
                     0.222289900];
    basis(2).coef = [-0.09996723
                     0.39951283
                     0.70011547];
    basis(3).name = 'P';
    basis(3).expo = [2.941249400
                     0.683483100
                     0.222289900];
    basis(3).coef = [0.15591627
                     0.60768372
                     0.39195739];                   
  case 'N'
    basis(1).name = 'S';
    basis(1).expo = [99.106169000
                     18.052312000
                     4.885660200];
    basis(1).coef = [0.15432897
                     0.53532814
                     0.44463454];
    basis(2).name = 'S';
    basis(2).expo = [3.780455900
                     0.878496600
                     0.285714400];
    basis(2).coef = [-0.09996723
                     0.39951283
                     0.70011547];
    basis(3).name = 'P';
    basis(3).expo = [3.780455900
                     0.878496600
                     0.285714400];
    basis(3).coef = [0.15591627
                     0.60768372
                     0.39195739];                   
  case 'O'
    basis(1).name = 'S';
    basis(1).expo = [130.7093200
                     23.8088610
                     6.4436083];
    basis(1).coef = [0.15432897
                     0.53532814
                     0.44463454];                
    basis(2).name = 'S';
    basis(2).expo = [5.0331513 
                     1.1695961 
                     0.3803890];
    basis(2).coef = [-0.09996723
                     0.39951283
                     0.70011547];
    basis(3).name = 'P';
    basis(3).expo = [5.0331513
                     1.1695961
                     0.3803890];
    basis(3).coef = [0.15591627
                     0.60768372
                     0.39195739];
  otherwise
    disp('STO-3G basis is not implemented for this atom');
end


end



