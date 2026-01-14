function [Gb,dR,DOFr,DOFm,AREA,R,DOFA,DOFB] = ...
    PERIODIC_BEAMS_MASS_MATRIXfun(a_A,a_B,DOMAINVAR,COOR,CONNECTb,TypeElementB,...
    DATA)

if nargin == 0 
    load('tmp1.mat')
end


ndim = size(COOR,2); 
da = a_A-a_B ;
%%%% FACE 1
iface=1 ; 
nodesfA = DOMAINVAR.NODES_faces12{1,iface} ; 
if  max(nodesfA-sort(nodesfA)) ~=0 
    error(['nodesfA should be sorted in ascending order'])
end
DOFA = small2large(nodesfA,ndim) ;   % DOFS face 1
COOR_FACE = COOR(nodesfA,:) ; % Coordinates of this face
[CentroidFA,AREA,Mst] =CentroidGeometricMassMatrixNEW(COOR,nodesfA,CONNECTb{iface},TypeElementB) ; 
COORrelA = bsxfun(@minus,COOR_FACE',CentroidFA')'; % Coordinates relative to centroid
BasisRrb = ConstructBasisRigidBody(COORrelA) ;
%% Mass interface matrix 
M = sparse(size(Mst,1)*ndim,size(Mst,1)*ndim); 
for idim =1:ndim
    INDLOC =idim:ndim:size(BasisRrb,1) ; 
    M(INDLOC,INDLOC) =  Mst ; 
end

if length(a_A) ==7
    % We have to add the warping modes
    Rwarping = ConstructWarpingMode(COORrelA) ;
    DATA = DefaultField(DATA,'MakeWarpingMode_Morthogonal',0) ; 
    if DATA.MakeWarpingMode_Morthogonal == 1
        % Make it M-orthogonal
        % In principle, it is not necessary (the mode it is by construction orthogonal)
    coeff = (BasisRrb'*M*BasisRrb)\(BasisRrb'*M*Rwarping) ;
    Rwarping = Rwarping - BasisRrb*coeff ;
    
    end
    R = [BasisRrb,Rwarping] ;
else
    R = BasisRrb; 
end

 
% FACE 2
iface=2 ; 
nodesfB = DOMAINVAR.NODES_faces12{end,iface} ;   % They are already paired,  
DOFB = small2large(nodesfB,ndim) ; % Set of DOFs face B
% Rorth = R'*M --->  R(INDX)'*N'*W*N 
Rbar = M*R ; 

% = zeros(size(R')) ; 
% for idim =1:ndim
%     INDLOC =idim:ndim:size(R,1) ; 
%     Rorth(:,INDLOC) = R(INDLOC,:)'*Mst ; 
% end

[Gb,dR,DOFr,DOFm] = BCs_BEAMS_PERIODIC(DOFA,DOFB,R,a_A,a_B,Rbar) ; 
 

 