function [Nshape,COORnodes,OTHER_OUTPUT] =  LinearLinearSolid_EIFEM(INFO_INTERFACE_MESH,CENTROID,COORbnd)
% JAHO, 6-APRIL-2025, SUNDAY, BALMES 185, BARCELONA
% Computation of Nshape for the interface ficti. modes, in EIFEM
% mIXED QUADRATIC/LINEAR INTERPOLATIONS, 2D. SOLID
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/107_EIFEM_plast3D/08_GENERAL_MODES.mlx
% --------------------------------------------------------------------------------------
if nargin == 0
    load('tmp.mat')
end

OTHER_OUTPUT = [] ; 
COORnodes = INFO_INTERFACE_MESH.COOR; %(:,2:end) ;

  
nmodes = size(COORnodes,1) ; % Number of modes
npoints = size(COORbnd,1) ;

Nshape = zeros(npoints,nmodes) ;

% cOORDINATES REFERRED TO THE CENTROID
ndim = size(CENTROID,2);
for idim = 1:ndim
    COORnodes(:,idim) = COORnodes(:,idim) - CENTROID(idim) ;
end
 


% NORMALIZATION COORDINATES 
% -------------------------
Lx = norm(COORnodes(1,:)-COORnodes(2,:)) ; 
Ly = norm(COORnodes(2,:)-COORnodes(3,:)) ; 

xi = 2*COORbnd(:,1)/Lx;
eta = 2*COORbnd(:,2)/Ly;
 

L1_xi = 0.5 *(1- xi);
L2_xi = 0.5 *(xi + 1);

% Lagrange polynomials in eta (linear)
M1_eta = 0.5*(1 - eta);
M2_eta = 0.5*(1 + eta);

Nshape = zeros(length(xi),4);  

inode = 1; 
Nshape(:,inode) = L1_xi.*M1_eta ; 
inode = 2; 
Nshape(:,inode) = L2_xi.*M1_eta ; 
inode = 3; 
Nshape(:,inode) = L2_xi.*M2_eta ; 
inode = 4; 
Nshape(:,inode) = L1_xi.*M2_eta ; 
 



% 
% for ielem = 1:length(ElemBnd)
%     
%     if length(ElemBnd{ielem}) == 2
%        % Linear interpolation         
%     Nshape =   LinearShapeFun2nodes(ElemBnd{ielem},COORnodes,COORbnd,Nshape) ; 
%         
%     elseif   length(ElemBnd{ielem}) == 3
%         % Quadratic interpolation
%       Nshape =   QuadraticShapeFun3nodes(ElemBnd{ielem},COORnodes,COORbnd,Nshape) ;   
%         
%         
%     end
%     
%     
% end

end
