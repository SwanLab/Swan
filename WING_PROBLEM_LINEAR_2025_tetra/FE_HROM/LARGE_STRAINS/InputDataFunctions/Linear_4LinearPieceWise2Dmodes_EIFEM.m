function [Nshape,COORnodes,OTHER_OUTPUT] =  Linear_4LinearPieceWise2Dmodes_EIFEM(INFO_INTERFACE_MESH,CENTROID,COORbnd)
% JAHO, 19-APRIL-2025, SATURDAY, Molinos Marfagoes, Cartagena
% Computation of Nshape for the interface ficti. modes, in EIFEM
%  LINEAR INTERPOLATIONS, 2-nodes two opposite edges, in 2D/ 4-nodes other edges.
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/107_EIFEM_plast3D/08_GENERAL_MODES.mlx
% --------------------------------------------------------------------------------------
if nargin == 0
    load('tmp.mat')
end
OTHER_OUTPUT = [] ;  


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




for ielem = 1:length(ElemBnd)
    
    if length(ElemBnd{ielem}) == 2
       % Linear interpolation         
    Nshape =   LinearShapeFun2nodes(ElemBnd{ielem},COORnodes,COORbnd,Nshape) ; 
        
    elseif   length(ElemBnd{ielem}) == 4
        % Quadratic interpolation
      Nshape =   LinearShapeFun4nodes(ElemBnd{ielem},COORnodes,COORbnd,Nshape) ;   
        
        
    end
    
    
end

end
