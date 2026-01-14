function[ Belem, wSTs, wST, XeALL, posgp,weigREP] = ComputeBelemALL(COOR,CN,TypeElement,nstrain,DATA) 
%--------------------------------------------------------------------------
% FUNCTION: ComputeBelemALL
%
% PURPOSE:
%   Computes the elemental strain-displacement matrices (`Belem`), as well as
%   associated Gauss integration data required in finite element formulations.
%
%   Specifically, it returns:
%     1) `Belem`: The global assembled B-matrix for all elements and Gauss points.
%                 This matrix maps nodal displacements to strains (or deformation gradients).
%     2) `wSTs` : Vector containing the product of Gauss weights and Jacobian determinants
%                at each Gauss point (used for integration).
%     3) `wST`  : Vector `wSTs` expanded to match the number of strain components.
%     4) `XeALL`: Coordinates of all elements, organized for vectorized operations.
%     5) `posgp`: Gauss point positions in the reference element.
%     6) `weigREP`: Replicated weights for each element (used in integration).
%
% INPUTS:
%   - COOR       : (nnode × ndim) Global coordinate matrix.
%   - CN         : (nelem × nnodeE) Connectivity matrix.
%   - TypeElement: (string) Finite element type (e.g., 'Quadrilateral').
%   - nstrain    : Number of strain components (e.g., 3 in 2D, 6 in 3D).
%   - DATA       : (struct) Additional settings, may include:
%       • DATA.DEFORMATION_GRADIENT == 1 → Return operator for deformation gradient (not strain).
%       • DATA.StrainStressWith4Components == 1 → Include ε_zz in 2D plane strain cases.
%       • DATA.posgp_given, DATA.weights_given → User-defined integration rule.
%
% OUTPUTS:
%   - Belem      : (nstrain × nelem × ngaus) × (ndim × nnodeE) Global B-matrix.
%   - wSTs       : (nelem × ngaus) Weights × det(Jacobian) per Gauss point.
%   - wST        : (nstrain × nelem × ngaus) Like `wSTs`, but expanded per component.
%   - XeALL      : (nelem × ndim) × nnodeE Element coordinate blocks.
%   - posgp      : (ndim × ngaus) Gauss point positions in reference element.
%   - weigREP    : (nelem × ngaus) Replicated weights.
%
% ALGORITHM:
%   - For each Gauss point:
%     1. Compute shape function derivatives in the parent domain.
%     2. Map derivatives to the physical domain using the Jacobian.
%     3. Transform to symmetric gradient (B-matrix for strain) or deformation gradient form.
%     4. Assign rows to the global `Belem` matrix.
%     5. Compute integration weights (`wST`, `wSTs`).
%
% SPECIAL OPTIONS:
%   - If `DATA.DEFORMATION_GRADIENT == 1`, returns 9 (3D) or 4 (2D) components of F.
%   - If `StrainStressWith4Components == 1`, adds dummy ε_zz component in 2D plane strain.
%
% NOTES:
%   - Highly vectorized for performance.
%   - This function is central to both small-strain and large-strain finite element formulations.
%   - It is compatible with both standard and custom quadrature rules.
%
% AUTHOR:
%   Joaquín A. Hernández Ortega
%   jhortega@cimne.upc.edu
%   CIMNE – Universitat Politècnica de Catalunya
%   First release: 26-Oct-2015
%   Revisions: Nov-2020, Feb-2022 (plane strain, deformation gradient)
%   Comments by ChatGPT4, 13-May-2025
%--------------------------------------------------------------------------

%%%%
% This subroutine   returns
%  1) the matrix Belem (nstrain*nelem*ngaus x ndime*nnodeE)  consisting of all element B-matrices
%  2) wSTs = zeros(nelem*ngaus,1) ;  % Vector containinig the product of weights and Jacobians at all gauss points
%  3) wST = zeros(nelem*ngaus*nstrain,1) ; % Similar to the WeightST, but replicating   each entry nstrain times
%  4)XeALL:  COORDINATES of all elements, arrangedin a nelem*ndim x nnodeE matrix
% Inputs:   COOR: Coordinate matrix (nnode x ndim), % CN: Connectivity matrix (nelem x nnodeE),
% TypeElement: Type of finite element (quadrilateral,...),
%% Vectorized version
% Joaquín A. Hernández (jhortega@cimne.upc.edu), 26-Oct-2015
%dbstop('10')
if nargin == 0
    load('tmp.mat')
elseif nargin == 4
    DATA = [] ;
end
nnode = size(COOR,1); ndim = size(COOR,2); nelem = size(CN,1); nnodeE = size(CN,2) ;
% nstrain = size(celasglo,1) ;
% Shape function routines (for calculating shape functions and derivatives)
DATA = DefaultField(DATA,'posgp_given',[]) ; 
if isempty(DATA.posgp_given)
TypeIntegrand = 'K';
else
   TypeIntegrand = {DATA.posgp_given,DATA.weights_given} ;  
end
%dbstop('20')
[weig,posgp,shapef,dershapef] = ComputeElementShapeFun(TypeElement,nnodeE,TypeIntegrand) ;
ngaus = length(weig) ;
% Initialization

% 20-Nov-2020 
%-------------
% COMPUTATION OF DEFORMATION GRADIENT VECTOR 
% ----------------------------------------------
% See /home/joaquin/Desktop/CURRENT_TASKS/POTENTIAL_RESEARCH_TOPICS/COMBINING_MULTISCALE_REDUCTIONMODELS/REPORT_MULTIS_REDUC_MODEL/
% NEW_IDEAS/DOCS/RAUL_BRAVO_THESHIS/ROTATION_SVD_20_Nov_2020/Rotation_SVD.pdf,
% page 33 onwards 
DATA = DefaultField(DATA,'DEFORMATION_GRADIENT',0) ; 
if DATA.DEFORMATION_GRADIENT==1
    % Solid elements 
    if nstrain == 3  ||  nstrain == 4
        % nstrain = 4, added 8-Feb-2022. Contemplate the case of plane
        % strain J2 plasticy
        % Fv1 = F11, Fv2 = F22, Fv3 =F12, Fv4 = F21. 
        % Fv = identity + Be*de 
        nstrain = 4 ; 
    elseif nstrain == 6
        % Fv1 = F11; Fv2 = F22; Fv3 = F33 ; Fv4 = F23  ; Fv5 = F13 ; Fv6 = F12  ; Fv7 = F32  ;
        % Fv8 = F31  ; Fv9  = F21 
        nstrain = 9; 
    else
        error('Option not implemented yet')
    end
       
else
    
end





Belem = zeros(nstrain*nelem*ngaus,nnodeE*ndim) ;
wSTs = zeros(nelem*ngaus,1) ;  % Vector containinig the product of weights and Jacobians at all gauss points
wST = zeros(nelem*ngaus*nstrain,1) ; % Similar to the WeightST, but replicating the each entry nstrain times
% COORDINATE MATRIX arranged in a nelem*ndim x nnodeE matrix
XeALL= COORallELEM(ndim,nelem,nnodeE,CN,COOR) ;
% Let us define a matrix ROWSgauss such that   Belem(ROWSgauss(g,:),:) returns
% the B-matrices of the g-th points of all elements
indREF = 1:nstrain*ngaus*nelem ;
ROWSgauss = reshape(indREF,nstrain,nelem*ngaus) ;
weigREP = repmat(weig',nelem,1)  ; % nelem x 1 tiling copies of weig

%DATA = DefaultField(DATA,'ISPSTRAIN',0) ;
DATA = DefaultField(DATA,'StrainStressWith4Components',0) ;

if  DATA.StrainStressWith4Components == 1 && nstrain ~=6  % && DATA.ISPSTRAIN ==1  
    % Plane strain including the nonzero component of strains (strain_z)
    INCLUDE_eZ = 1;
else
    INCLUDE_eZ = 0;
end



for  g = 1:ngaus
    % Matrix of derivatives for Gauss point "g"
    BeXi = dershapef(:,:,g) ;
    % Jacobian Matrix for the g-th G. point of all elements %
    JeALL = XeALL*BeXi' ;
    %%%%%%%%%
    % JAcobian
    detJeALL= determinantVECTORIZE(JeALL,ndim) ;
    % Matrix of derivatives with respect to physical coordinates
    inv_JeTall = inverseTRANSvectorize(JeALL,ndim,detJeALL) ;
    BeTILDEall = inv_JeTall*BeXi ;
    % Matrix of symmetric gradient
    % dbstop('18')
    if DATA.DEFORMATION_GRADIENT == 1
         BeALL = QtransfBvect_F(BeTILDEall,ndim) ;
    else
        if INCLUDE_eZ == 0
            BeALL = QtransfBvect(BeTILDEall,ndim) ; % Transformation from B-matrix for scalar fields to B-matrix for vector fields
        else
            BeALL = QtransfBvect_2Dps(BeTILDEall,ndim) ; % Plane strain, including zeros for e_z
        end
    end
    
    % Assigning BeALL ( g-th Gauss point) to the global matrix Belem
    ROWSglo =  ROWSgauss(:,g:ngaus:ngaus*nelem);
    ROWSglo = ROWSglo(:) ;
    Belem(ROWSglo,:) = BeALL ;
    % Weight vectors
    % --------------
    [wST,wSTs] = DetermineWeightsST(detJeALL,weigREP,ngaus,g,nelem,wST,wSTs,nstrain,ROWSglo);
    
    
    
end
