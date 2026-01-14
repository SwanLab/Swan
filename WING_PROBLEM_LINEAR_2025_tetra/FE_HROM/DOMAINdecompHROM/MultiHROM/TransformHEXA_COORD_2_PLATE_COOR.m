function Vplate = TransformHEXA_COORD_2_PLATE_COOR(Vfe,DeltaZ) ; 
% Given the 8x3 = 24 modes associated to the shape functions "modes" of a
% cube (Vfe), this function returns the constrained modes arising from
% imposing the "plate" constraints (u_z =0). Vplate is a matrix of 20
% columns
% Vplate  [Vplate_corner1, Vplate_corner2,Vplate_corner3,Vplate_corner4]
% See 
%  JAHO 16-Feb-2023
% ----------------------------------------------------------------------
if nargin == 0
    load('tmp1.mat')
end
ndofs =size(Vfe,1) ; 
nmodes_constr = 20 ; 
ndim = 3; 
 nnode = 4;  
 % INDEXES NODES HEXAHEDRA
 % ---------------------------------------

zM = ndim:ndim:ndim*nnode; % INDEXES Z-DISPLACEMENTS BOTTOM SURFACE 
zP = ndim*nnode + (ndim:ndim:ndim*nnode); % INDEXES Z-DISPLACEMENTS top SURFACE 
xM = 1:ndim:ndim*nnode; % INDEXES x-DISPLACEMENTS BOTTOM SURFACE 
xP = ndim*nnode + (1:ndim:ndim*nnode); % INDEXES x-DISPLACEMENTS top SURFACE
yM = 2:ndim:ndim*nnode; % INDEXES x-DISPLACEMENTS BOTTOM SURFACE 
yP = ndim*nnode + (2:ndim:ndim*nnode); % INDEXES x-DISPLACEMENTS top SURFACE 
Vzm = Vfe(:,zM) + Vfe(:,zP) ;  
Vplate = zeros(ndofs,nmodes_constr) ; 

ndimP = 5; 
% Translations in the z-direction 
z = 3:ndimP:nmodes_constr ; 
Vplate(:,z) = Vzm ; 
for inode= 1:nnode 
    indPLATE  = ((inode-1)*ndimP)+[1,2,4,5] ; 
      
    
    indHEXA  = [xM(inode),xP(inode),yM(inode),yP(inode)] ; 
    B = [0.5 0.5 0 0 
        0  0  0.5  0.5 
        0  0   1/DeltaZ(inode)  -1/DeltaZ(inode)
        -1/DeltaZ(inode)  1/DeltaZ(inode) 0 0 ] ;  
        
    A = inv(B) ;     
    
    Vplate(:,indPLATE)  = Vfe(:,indHEXA)*A ; 
end