function [Gb,dR,DOFr,DOFm,AREA,R,DOFA,DOFB] = ...
    PRESCRIBED_FLUCTUATIONS_withFORCESonly(U_b,DOMAINVAR,COOR,CONNECTb,TypeElementB,...
    DATA,gammaB)
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/105_BackPeriodic/03_IMPOSED_FORCES.mlx
% JAHO, 9-Nov-2023, Barcelona, UPC, CIMNE
if nargin == 0 
    load('tmp.mat')
end


ndim = size(COOR,2); % Number of spatial dimensions
%%%% FACE A 
iface=1 ; 
nodesfA = DOMAINVAR.NODES_faces12{1,iface} ; % Node indexes
if  max(nodesfA-sort(nodesfA)) ~=0 
    error(['nodesfA should be sorted in ascending order'])
end
DOFA = small2large(nodesfA,ndim) ;   % DOFS face 1
COOR_FACE = COOR(nodesfA,:) ; % Coordinates of this face
[CentroidFA,AREA,Mst] =CentroidGeometricMassMatrixNEW(COOR,nodesfA,CONNECTb{iface},TypeElementB) ; 
COORrelA = bsxfun(@minus,COOR_FACE',CentroidFA')'; % Coordinates relative to centroid
BasisRrb = ConstructBasisRigidBody(COORrelA) ; % Rigid body modes 
%% Mass interface matrix 
M = sparse(size(Mst,1)*ndim,size(Mst,1)*ndim); 
for idim =1:ndim
    INDLOC =idim:ndim:size(BasisRrb,1) ; 
    M(INDLOC,INDLOC) =  Mst ; 
end

% if length(a_A) ==7
%     % We have to add the warping modes
%     Rwarping = ConstructWarpingMode(COORrelA) ;
%     DATA = DefaultField(DATA,'MakeWarpingMode_Morthogonal',0) ; 
%     if DATA.MakeWarpingMode_Morthogonal == 1
%         % Make it M-orthogonal
%         % In principle, it is not necessary (the mode it is by construction orthogonal)
%     coeff = (BasisRrb'*M*BasisRrb)\(BasisRrb'*M*Rwarping) ;
%     Rwarping = Rwarping - BasisRrb*coeff ;
%     
%     end
%     R = [BasisRrb,Rwarping] ;
% else
     R = BasisRrb; 
% end

 
% FACE 2
iface=2 ; 
nodesfB = DOMAINVAR.NODES_faces12{end,iface} ;   % They are already paired,  
DOFB = small2large(nodesfB,ndim) ; % Set of DOFs face B
% Rorth = R'*M --->  R(INDX)'*N'*W*N 
 Rbar = M*R ; 
 % 
% if  ~isempty(U_b)
    
[Gb,dR,DOFr,DOFm] = BCs_PRESCRIBED_FLUCTUATIONS_withFORCESonly(DOFA,DOFB,R,U_b,gammaB,Rbar,M) ; 
% else
%    
%     [Gb,dR,DOFr,DOFm] = BCs_PRESCRIBED_FLUCTUATIONS_withFORCES_nobend(DOFA,DOFB,R,U_n) ; 
% 
% end
 

 