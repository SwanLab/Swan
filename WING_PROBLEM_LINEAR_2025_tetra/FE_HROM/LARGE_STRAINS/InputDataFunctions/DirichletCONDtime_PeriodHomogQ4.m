function [DOFr,dR,DISP_CONDITIONS,INFO_PERIODIC_CONDITIONS] = DirichletCONDtime_PeriodHomogQ4(DIRICHLET,DATA,ndim,MESH,GEOproperties,DATALOC)
%
%  
%
% JAHO- 7-Oct-2O25, Upc Campus TErrassa 
%--------------------------------------------------------------------------------------
if nargin == 0
    load('tmp.mat')
end
INFO_PERIODIC_CONDITIONS = [] ;

% The following are linear boundary conditions
[DOFr_LINEAR,dR_LINEAR] = DirichletCONDtime_homogZERO(DIRICHLET,DATA,ndim,MESH,GEOproperties,DATALOC) ; 
% How to adapt now to PERIODIC boundary conditions ? 
 
DATALOC  = DefaultField(DATALOC,'LabelEntitiesDefiningBoundary',[1,2,3,4]) ;
 FACES_BND  = DATALOC.LabelEntitiesDefiningBoundary ; % FACES DEFINING THE BOUNDARIES OF THE FINE MESH DOMAIN
 NODES_FACES = MESH.NODES_FACES(FACES_BND) ;  %
 NODES_FACES = NODES_FACES(:) ;               %
 rnodLOC = unique(cell2mat(NODES_FACES)) ;  % LIST OF NODES WITH PRESCRIBED DISPLACEMENTS (BOUNDARY)
% % ----------------------------------------------------------
 COORbnd = MESH.COOR(rnodLOC,:) ;  % COORDINATES
 
  
[ DISP_CONDITIONS,INFO_PERIODIC_CONDITIONS] = PeriodicQ4_homog(DIRICHLET,DATA,ndim,MESH,GEOproperties,DATALOC,DOFr_LINEAR,dR_LINEAR);
        
dR=  [] ; 
DOFr = [] ; 