function [Nshape,COORnodes,OTHER_OUTPUT] =  LinearQuadPieceWise2Dmodes_EIFEM(INFO_INTERFACE_MESH,CENTROID,COORbnd)
% JAHO, 5-APRIL-2025, SATURDAY, BALMES 185, BARCELONA
% Computation of Nshape for the interface ficti. modes, in EIFEM
% mIXED QUADRATIC/LINEAR INTERPOLATIONS, 2D.
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/107_EIFEM_plast3D/08_GENERAL_MODES.mlx
% --------------------------------------------------------------------------------------
if nargin == 0
    load('tmp.mat')
end


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
    elseif   length(ElemBnd{ielem}) == 3
        % Quadratic interpolation
         [Nshape,xiEDGES{ielem},IndPointsBNDedge{ielem}]  =   QuadraticShapeFun3nodes(ElemBnd{ielem},COORnodes,COORbnd,Nshape) ;
        
         NumberNodesPerEdge(ielem) = 3;
    end
    
    
end

OTHER_OUTPUT.xiEDGES = xiEDGES ;
OTHER_OUTPUT.IndPointsBNDedge = IndPointsBNDedge ;
OTHER_OUTPUT.NumberNodesPerEdge = NumberNodesPerEdge ;
