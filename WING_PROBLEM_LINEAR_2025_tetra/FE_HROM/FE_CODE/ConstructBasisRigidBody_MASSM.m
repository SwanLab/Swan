function [BasisRrb,VOLUME,Rbar,Centroid,INERTIA ]= ConstructBasisRigidBody_MASSM(COORabs,Nst,wSTs,DATAIN)
%--------------------------------------------------------------------------
% FUNCTION: ConstructBasisRigidBody_MASSM
%
% PURPOSE:
%   Computes the **rigid body modes**, **volume (or area)**, **centroid**, and
%   **rotational inertia tensor** of a finite element subdomain, using mass-like
%   integration over the domain via the shape functions `Nst` and weights `wSTs`.
%   This function is commonly used in the context of reduced-order modeling (ROM)
%   and multiscale methods (e.g., EIFEM), where rigid body modes form the basis
%   for describing non-deformational displacements.
%
% INPUTS:
%   - COORabs : (nnode × ndim) matrix of absolute coordinates of the domain's nodes.
%   - Nst     : (ngauss_total × nDOF) shape function matrix at Gauss points.
%   - wSTs    : (ngauss_total × 1) vector of integration weights times Jacobian determinants.
%   - DATAIN  : Structure with optional field:
%       * .ORTHOGONAL_RIGID_BODY_MODES :
%           0 → no orthogonalization (default)
%           1 → orthogonalize using QR decomposition
%           2 → normalize each mode individually
%
% OUTPUTS:
%   - BasisRrb : (nDOF × nRBM) matrix of rigid body modes (columns: 3 in 2D, 6 in 3D)
%   - VOLUME   : Scalar representing the volume (or area in 2D) of the domain.
%   - Rbar     : Matrix = Nstᵗ * diag(wSTs) * (Nst * BasisRrb), used in projection.
%   - Centroid : (1 × ndim) vector, geometric centroid of the domain.
%   - INERTIA  : Rotational inertia matrix of the domain, extracted from Rbar.
%                * 3×3 in 3D
%                * scalar in 2D (entry (3,3) of GEOMETRIC_PROPERTIES_VOLUME)
%
% DETAILS:
%   - Rigid body modes are constructed analytically based on the relative
%     coordinates of each node w.r.t. the centroid.
%       In 2D: 2 translations + 1 rotation
%       In 3D: 3 translations + 3 rotations (using antisymmetric spin matrices)
%
%   - The centroid is computed as the mass-weighted average of nodal coordinates:
%         Centroid = (∫ N(x) x dx) / (∫ dx)
%
%   - Inertia tensor is assembled from the projected mass matrix of the rigid modes.
%
% APPLICATION:
%   - Construction of rigid body basis functions for domain decomposition methods
%   - Enforcing kinematic admissibility in multi-subdomain ROM
%   - Filtering out rigid body motion in stress computations
%
% SEE ALSO:
%   GeometricVarDOMAINS, ComputeBelemALL, AssemblyNGlobal, QtransfBvect
%
% AUTHOR:
%   Joaquín A. Hernández Ortega (CIMNE/UPC)
%   Last updated: 28-Jan-2025, Barcelona
%   Comments by ChatGPT4, 13-May-2025
%--------------------------------------------------------------------------

if nargin == 0
    load('tmp.mat')
end


%%%%% COMPUTING THE CENTROID 


ndim = size(COORabs,2) ; 
wSTdiag = CompWeightDiag(wSTs,ndim)  ;
%wSTdim = repmat(wSTs',ndim,1)' ; 
%Mst = (wSTdiag*Nst)'*Nst ;
% Recomputing centroid
Centroid = zeros(1,ndim) ;
VOLUME = sum(wSTs) ;
%for idim = 1:3
COOR = COORabs' ; 
NstCOOR = Nst*COOR(:) ; 
Centroid = zeros(1,ndim) ; 
for idim = 1:ndim
    Centroid(idim) = wSTs'*NstCOOR(idim:ndim:end)/VOLUME ;
end
% 
COORrel = bsxfun(@minus,COORabs',Centroid')'; % Coordinates relative to centroid

if ndim == 2
    BasisRrb =  zeros(prod(size(COORrel)),3) ;
    BasisRrb(1:ndim:end,1) = 1;
    BasisRrb(2:ndim:end,2) = 1;
    BasisRrb(1:ndim:end,3) = -COORrel(:,2);
    BasisRrb(2:ndim:end,3) = +COORrel(:,1);  % or the other way around ?
else
    
    BasisRrb =  zeros(prod(size(COORrel)),6) ;
    BasisRrb(1:3:end,1) = 1 ;
    BasisRrb(2:3:end,2) = 1;
    BasisRrb(3:3:end,3) = 1;
    
  
    
      % ------------------------------------------------
    % CHANGE INTRODUDUCED 3th-october-2021
    %------------------------------------
    % OLD VERSION 
    % ---------------------------
%      % Rotation modes
  % OLD VERSION
%     BasisRrb(2:3:end,4) = COORrel(:,3);
%     BasisRrb(3:3:end,4) = -COORrel(:,2);
%     
%     BasisRrb(1:3:end,5) = -COORrel(:,3);
%     BasisRrb(3:3:end,5) = COORrel(:,1);   
%     
%       
%     BasisRrb(1:3:end,6) = COORrel(:,2);
%     BasisRrb(2:3:end,6) = -COORrel(:,1);
%   

% NEW VERSION 
% ----------------------------
% See  /home/joaquin/Desktop/CURRENT_TASKS/PAPERS_2020_onwards/MultiECMgen/DiscussionRIGIDBODY.xoj  /*.tex 
    % Rotation modes
    % [R_r]_I   = ROtation modes of node I, of coordinates [X]_I = [x1,x2,x3]
    %                MODE  4   5   6
    %                ------------------
    % = -spin([X]_I)  =  [0  +x3 -x2 
    %                      -x3  0  +x1
    %                      +x2  -x1  0] = 
    
    % Rotation modes
  
     BasisRrb(2:3:end,4) = -COORrel(:,3);
    BasisRrb(3:3:end,4) = +COORrel(:,2);
    
    BasisRrb(1:3:end,5) = +COORrel(:,3);
    BasisRrb(3:3:end,5) = -COORrel(:,1);
 
     BasisRrb(1:3:end,6) = -COORrel(:,2);
     BasisRrb(2:3:end,6) = +COORrel(:,1);
   
    
    
    
    
end

DATAIN = DefaultField(DATAIN,'ORTHOGONAL_RIGID_BODY_MODES',0) ;

if DATAIN.ORTHOGONAL_RIGID_BODY_MODES== 1
    BasisRrb = orth(BasisRrb)  ;
elseif DATAIN.ORTHOGONAL_RIGID_BODY_MODES== 2
    NNN = sqrt(sum(BasisRrb.^2,1)) ;
    for imode = 1:size(BasisRrb,2)
        BasisRrb(:,imode) = BasisRrb(:,imode)/NNN(imode) ;
    end
end

Rbar = Nst'*wSTdiag*(Nst*BasisRrb) ; 

GEOMETRIC_PROPERTIES_VOLUME = (BasisRrb'*Rbar) ;  % Volume and moment of inertias

if ndim == 3
    INERTIA =GEOMETRIC_PROPERTIES_VOLUME(4:6,4:6) ;
  %  MSG{end+1} =['VOLUME =',num2str(I(1,1))] ;
   % MSG{end+1} ='MATRIX OF INERTIAS ';
   % MSG{end+1} ='------------------------------' ;
   % MSG{end+1} =num2str(Inert,3)  ;
else
   % MSG{end+1} =['Area (slice) = ',num2str(I(1,1))]; ;
   % MSG{end+1} =['Moment of Inertia (slice) = ',num2str(I(3,3))]  ;
    
    INERTIA = GEOMETRIC_PROPERTIES_VOLUME(3,3) ;
end


%PROPERTIES = Rbar'*BasisRrb;   This matrix contains 