function [Nshape,COORnodes,OTHER_OUTPUT] =  LinearQuarticPieceWise2Dmodes_EIFEM(INFO_INTERFACE_MESH,CENTROID,COORbnd)
% JAHO, 8-MAY-2025, THURSDAY, upc TERRASSA
% Computation of Nshape for the interface ficti. modes, in EIFEM
% mIXED QUARTIC/LINEAR INTERPOLATIONS, 2D.
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/12_REPETITIVE_TRAIN.mlx
% --------------------------------------------------------------------------------------
if nargin == 0
    load('tmp.mat')
end
OTHER_OUTPUT = [] ;

INFO_INTERFACE_MESH = DefaultField(INFO_INTERFACE_MESH,'xi_nodes',[]) ; 
xi_nodes = INFO_INTERFACE_MESH.xi_nodes; 

COORnodes = INFO_INTERFACE_MESH.COOR; %(:,2:end) ;
ElemBnd= INFO_INTERFACE_MESH.LINES ;

nmodes = size(COORnodes,1) ; % Number of modes
npoints = size(COORbnd,1) ;

Nshape = zeros(npoints,nmodes) ;

% cOORDINATES REFERRED TO THE CENTROID
ndim = size(CENTROID,2);
for idim = 1:ndim
    COORnodes(:,idim) = COORnodes(:,idim) - CENTROID(idim) ;
end

xiEDGES = cell(size(ElemBnd)) ;
IndPointsBNDedge = cell(size(ElemBnd)) ;
NumberNodesPerEdge = zeros(size(ElemBnd)) ;


for ielem = 1:length(ElemBnd)
    
    if length(ElemBnd{ielem}) == 2
        % Linear interpolation
        [Nshape,xiEDGES{ielem},IndPointsBNDedge{ielem}] =   LinearShapeFun2nodes(ElemBnd{ielem},COORnodes,COORbnd,Nshape) ;
        NumberNodesPerEdge(ielem) = 2;
        
    elseif   length(ElemBnd{ielem}) == 5
        % Quadratic interpolation
        [Nshape,xiEDGES{ielem},IndPointsBNDedge{ielem}] =   QuarticShapeFun5nodes(ElemBnd{ielem},COORnodes,COORbnd,Nshape,xi_nodes) ;
        
        NumberNodesPerEdge(ielem) = 5;
        
    end
    
    
end
OTHER_OUTPUT.xiEDGES = xiEDGES ;
OTHER_OUTPUT.IndPointsBNDedge = IndPointsBNDedge ;
OTHER_OUTPUT.NumberNodesPerEdge = NumberNodesPerEdge ;
