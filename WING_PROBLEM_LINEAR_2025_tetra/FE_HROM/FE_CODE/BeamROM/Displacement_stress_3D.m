function [DISP3D,DISP3D_lateral,stressGLO,STRESS_DATA,DATAIN,DATA_REFMESH] ...
    = Displacement_stress_3D(DATAROM,MESH1D,qDEF,qRB,DATAIN,DATA_REFMESH)


if nargin ==0
    load('tmp.mat')
end

% SELECTED_DOMAINS = [] ;
% if length(DATAROM) == 1
%     [DISP3D,DISP3D_lateral,stressGLO,STRESS_DATA,DATAIN,DATA_REFMESH] ...
%     = Displacement_stress_3D_1slice(DATAROM{1},MESH1D,qDEF,qRB,DATAIN,DATA_REFMESH{1}) ;
%
% else
     [DISP3D,DISP3D_lateral,stressGLO,STRESS_DATA,DATAIN,DATA_REFMESH] ...
    = Displacement_stress_3D_JOINTslice(DATAROM,MESH1D,qDEF,qRB,DATAIN,DATA_REFMESH) ;

%end
