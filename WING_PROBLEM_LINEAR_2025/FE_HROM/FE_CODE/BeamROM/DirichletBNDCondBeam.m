function [DOFr,DOFl,dR] = DirichletBNDCondBeam(DATAROM,MESH1D,DISP,DATAIN,ndim)

if nargin == 0
    load('tmp1.mat')
end

if length(DATAROM) == 1
    % Just one type of slice
    [DOFr,DOFl,dR] = DirichletBNDCondBeam_onlybeam(DATAROM{1},MESH1D,DISP,DATAIN,ndim) ; 
    
else
    % Slices and joints, but with just 6 modes
      [DOFr,DOFl,dR] = DirichletBNDCondBeam_mixed(DATAROM{1},MESH1D,DISP,DATAIN,ndim) ; 
    
end
 