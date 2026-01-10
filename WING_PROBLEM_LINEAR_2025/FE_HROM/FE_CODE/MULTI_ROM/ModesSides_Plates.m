function [Vside,As1,As2 ]= ModesSides_Plates(iSIDE,DATA_REFMESH)

if nargin == 0
    load('tmp.mat')
end


% %
% nodesSIDE= DATA_REFMESH.NODES_SIDES{iSIDE} ;
% COOR_FACE = DATA_REFMESH.COOR(nodesSIDE,:) ; % Coordinates of this SIDE
% CentroidFA =  DATA_REFMESH.CENTRf{iSIDE} ; % Centroid of the FACE (FACE = SIDE + CORNERS)
% COORrelA = bsxfun(@minus,COOR_FACE',CentroidFA')'; % Coordinates relative to centroid
% % Rigid body modes (6)
% R = ConstructBasisRigidBody(COORrelA) ;
 nodesFACE= DATA_REFMESH.NODES_faces{iSIDE} ;
 COOR_FACE = DATA_REFMESH.COOR(nodesFACE,:) ; % Coordinates of this SIDE
 CentroidFA =  DATA_REFMESH.CENTRf{iSIDE} ; % Centroid of the FACE (FACE = SIDE + CORNERS)
 COORrelA = bsxfun(@minus,COOR_FACE',CentroidFA')'; % Coordinates relative to centroid
 % Rigid body modes (6)
 R = ConstructBasisRigidBody(COORrelA) ;
% Deformational modes (4). We determine these modes by first computing the
% 10 modes  arising from a linear distribution
Vlinear= ModesLocSides(COORrelA,iSIDE) ;
% Now we inspect mode by mode to check which modes is linearly independent
% from the rigid body modes

V = R ;  
ModesIncluded = [] ;
for imode = 1:size(Vlinear,2)
    Vexp = [V Vlinear(:,imode)] ;
    if rank(Vexp) == size(Vexp,2)
        V = Vexp ;
        ModesIncluded = [ModesIncluded;imode] ;
    end
end
 
%
% Therefore, the desired basis is given by (rigid body + deformational)
Vface = V ;

 % if DECOMPOSITION == 0 
% Vface = V ; 
% 
% end

 
% This matrix includes both corner DOFs and side DOFs. Let us decompose it
% into corner and sides 
% Matrix of corner modes 
Vc = DATA_REFMESH.ModesCornes{1} ;  
nDOFScorner = size(Vc,1); % Number of DOFs associated to one corner
% Hence
Vc1c2 = Vface(1:2*nDOFScorner,:) ; 
% On the other hand, the desired matrix (side DOFs) is given  bu 
Vside = Vface((2*nDOFScorner+1):end,:) ;  

% Connection between the amplitudes associated to V and Vc 
VcDIAG = blkdiag(Vc,Vc) ; 
Aside_corners  = Vc1c2\VcDIAG ; 

As1 = Aside_corners(:,1:5) ; 
As2 = Aside_corners(:,6:end) ; 

% For checking the result 
Acorners_sides  = inv(Aside_corners) ; 



% C


 





end

function  V= ModesLocSides(COORrelA,iSIDE)


if iSIDE == 1
    iREF = 2;
elseif iSIDE == 2
    iREF = 1;
end


nmodes = 10 ;
ndim = size(COORrelA,2) ;

ndof = prod(size(COORrelA)) ;
V = zeros(ndof,nmodes) ;


% Face 1
% ---------

%% --------------------------
% x-entries
idim = 1;

imode = 1;
V(idim:ndim:end,imode) = 1 ;

imode = 2;
V(idim:ndim:end,imode) = COORrelA(:,iREF) ;

imode = 3;
V(idim:ndim:end,imode) = COORrelA(:,3) ;

imode = 4;
V(idim:ndim:end,imode) = COORrelA(:,3).*COORrelA(:,iREF) ;

% y-entries
idim = 2;

imode = 5;
V(idim:ndim:end,imode) = 1 ;

imode = 6;
V(idim:ndim:end,imode) = COORrelA(:,iREF) ;

imode = 7;
V(idim:ndim:end,imode) = COORrelA(:,3) ;

imode = 8;
V(idim:ndim:end,imode) = COORrelA(:,3).*COORrelA(:,iREF) ;

% z-entries
idim = 3 ;
imode = 9;
V(idim:ndim:end,imode) = 1 ;

imode = 10;
V(idim:ndim:end,imode) = COORrelA(:,iREF) ;




end