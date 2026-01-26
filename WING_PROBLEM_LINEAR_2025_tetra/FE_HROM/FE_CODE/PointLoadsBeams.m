function Fpnt =  PointLoadsBeams(DOMAINVAR,COOR,CONNECTb,...
    TypeElementB,FORCES,DATA)

if nargin == 0
    load('tmp2.mat')
end
% ---------
% FACE A
% ---------
iface = 1;
ndim = 3;
nodesfA = DOMAINVAR.NODES_faces12{1,iface} ;
DOFA = small2large(nodesfA,ndim) ;   % DOFS face 1

COOR_FACE = COOR(nodesfA,:) ; % Coordinates of this face
[CentroidFA,AREA,Mst] =CentroidGeometricMassMatrixNEW(COOR,nodesfA,CONNECTb{iface},TypeElementB) ;
COORrelA = bsxfun(@minus,COOR_FACE',CentroidFA')'; % Coordinates relative to centroid
R = ConstructBasisRigidBody(COORrelA) ;

M = sparse(size(Mst,1)*ndim,size(Mst,2)*ndim) ;
for idim = 1:ndim
    M(idim:ndim:end,idim:ndim:end) = Mst ;
end

b_A_input = FORCES{1}(:) ; 

Fpnt = zeros(size(COOR,1)*ndim,1) ;
b_A = (R'*M*R)\b_A_input ;
% Nodal forces FACE A
Fpnt(DOFA) = M*R*b_A ;
%
% FACE B
% ---------
iface = 2;
nodesfB = DOMAINVAR.NODES_faces12{end,iface} ;
DOFB = small2large(nodesfB,ndim) ;   % DOFS face 1

DATA = DefaultField(DATA,'angROTATION_FACE',[]) ;
DOMS = [1,size(CONNECTb,1)] ;

 if  ~isempty(DATA.angROTATION_FACE)
    total_angle = DOMS(iface)*DATA.angROTATION_FACE ;
    
     ROTMATRIX = LocalRotMatrix(iface,DATA,DOMS,ndim) ; 
    
%     ROTMATRIX = [cos(total_angle),sin(total_angle)
%         -sin(total_angle)  ,   cos(total_angle)] ;
%     if ndim == 3
%         ROTMATRIX_loc = eye(3) ;
%         ROTMATRIX_loc(1:2,1:2) = ROTMATRIX ;
%         ROTMATRIX = ROTMATRIX_loc ;
%     end
    %         uBARloc = ROTMATRIX*reshape(uBARloc,ndim,[]) ;
    %         uBARloc = uBARloc(:) ;
else
    ROTMATRIX = [];
    
 end
 

%COOR_FACE = COOR(nodesfB,:) ; % Coordinates of this face
%[CentroidFA,AREA,Mst] =CentroidGeometricMassMatrixNEW(COOR,nodesfA,CONNECTb{iface},TypeElementB) ;
%COORrelA = bsxfun(@minus,COOR_FACE',CentroidFA')'; % Coordinates relative to centroid

%   % ROTATION 
%     % -------
%     if ~isempty(ROTMATRIX)
%         COORrel = (ROTMATRIX'*COORrelA')' ; 
%     end
%     
% 
% R = ConstructBasisRigidBody(COORrel) ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -----------------------------------
%%%%% NODAL FORCES
% ----------------------------------

b_B_input = FORCES{2}(:) ; 


b_B = (R'*M*R)\b_B_input ;
floc = M*R*b_B  ; 
% Nodal forces FACE B
if   ~isempty(DATA.angROTATION_FACE)
    
    floc = ROTMATRIX*reshape(floc,ndim,[]) ;
    floc = floc(:) ;
    
end


Fpnt(DOFB) = floc ;
