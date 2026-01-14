function [stress,der_stress] = StressesFromStrainCABLE(MATPRO,e,DATA) 
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/FIBREGY_PROJECT_2022/04_MOORING_PROBLEMS/01_STATICgrav.mlx
if nargin == 0
    load('tmp1.mat')
end

stress = zeros(size(e)) ; 
der_stress = zeros(size(e)) ; 
if length(MATPRO.PROPMAT) > 1
    error('Adapt this part for cases with more than one material (for the hyperreduction problem )')
    
end

%for imat = 1:length(MATPRO.PROPMAT)    
 %   ELEMS = find(DATA.MESH1D.MaterialType == imat) ;
  %  DOFSelem = small2large(ELEMS,DATA.MESH.ngaus) ; 
  imat = 1; 
    [stress ,der_stress] = feval(MATPRO.PROPMAT(imat).stress_versus_strain,MATPRO.PROPMAT(imat).InputsCONSTEQ,e   ) ;
    
%     if DATA.CALC_CTANG
%         der_stress(DOFSelem)
%     end
    
%end