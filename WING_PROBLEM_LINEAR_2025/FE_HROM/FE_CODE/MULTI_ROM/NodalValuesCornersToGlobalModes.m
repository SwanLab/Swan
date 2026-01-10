function S = NodalValuesCornersToGlobalModes(NODES_CORNERS,x,y,z,COOR,Q)
% Matrix relating amplitude of corner modes (5 modes per corner) with
% global modes (trilinear displacement, hexagonal finite element)
% 
ndim = size(COOR,2); 

% -------------------------------------
% Vc: Matrix of modes of corner lines 
% ------------------------------------
CNODES = NODES_CORNERS{1} ;
% Centroid corner line = 1 
xC(1) = x.min ; 
xC(2) = y.max ; 
xC(3) = 0.5*(z.min+z.max) ; 
COOR_c = COOR(CNODES,:) ; 
COORrelA = bsxfun(@minus,COOR_c',xC')'; % Coordinates relative to centroid
R = ConstructBasisRigidBody(COORrelA) ;
Vc = R(:,1:5) ; 
% ----------------------------------------------
% ----------
% Solving system Qc*b = Vc*a  
% [Q(c1,:)',Q(c2,:)' ... ]*b = diag[Vc,Vc...]*a 
NDOFc = length(NODES_CORNERS{1})*ndim ; 
ncorner = 4; 
Qc = Q(1:(ncorner*NDOFc),:) ; 
% 
VcDIAG = cell(ncorner,1) ; 
VcDIAG(:) = {Vc} ; 
VcDIAG = blkdiag(VcDIAG{:}) ; 
% MAtrix b = S*a
S = Qc\VcDIAG ; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%