function [DOFr,dR,DISP_CONDITIONS,INFO_PERIODIC_CONDITIONS] = DirichletCONDtime_Period1D(DIRICHLET,DATA,ndim,MESH,DATALOC)
% Goal. Determine DOFr and    dR(t) .
% Periodicity conditions in one direction
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/108_EIFEM_metamat/01_meta_1D.mlx
%
% JAHO- 14-May-2O24, Terrassa,  UPC
%--------------------------------------------------------------------------------------
if nargin == 0
    load('tmp3.mat')
end
INFO_PERIODIC_CONDITIONS = [] ;

% % BOUNDARY NODES THAT WILL BE SUBJECT TO
% % THE FE-SHAPE-LIKE BOUNDARY CONDITIONS
% FACES_BND  = DATALOC.LabelEntitiesDefiningBoundary ; % FACES DEFINING THE BOUNDARIES OF THE FINE MESH DOMAIN
% NODES_FACES = MESH.NODES_FACES(FACES_BND) ;  %
% NODES_FACES = NODES_FACES(:) ;               %
% rnodLOC = unique(cell2mat(NODES_FACES)) ;  % LIST OF NODES WITH PRESCRIBED DISPLACEMENTS (BOUNDARY)
% % ----------------------------------------------------------
% COORbnd = MESH.COOR(rnodLOC,:) ;  % COORDINATES

 DOFr = [] ; dR= [] ; 
if ndim == 2    
    [ DISP_CONDITIONS,INFO_PERIODIC_CONDITIONS] = Periodic_1direct(DIRICHLET,DATA,ndim,MESH,DATALOC);    
    DOFr = DISP_CONDITIONS.DOFr ; 
    dR = DISP_CONDITIONS.dR ; 
else
    error('3D option not implemented yet ')
end



