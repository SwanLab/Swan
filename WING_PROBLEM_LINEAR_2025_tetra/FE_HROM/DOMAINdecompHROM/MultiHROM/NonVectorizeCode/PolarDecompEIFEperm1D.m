function [TRANS_COOR,normD] = PolarDecompEIFEperm1D(Xe,EIFEoper,DATA)
% Adaptation of  PolarDecompEIFEperm.m to 1D problems
% JAHO,13-May-2024
% ----------------------------------------------------------------------------
if nargin == 0
    load('tmp2.mat')
elseif  nargin == 2
    DATA = [] ;
    DATA.ToleranceDeformationalPartTransformationParentDomain = 1e-4 ;
end
[ndim,nnodeE] = size(Xe) ;

% COMPUTATION OF THE DILATION FACTOR (ACTUALLY THE DETERMINANT OF THE
% JACOBIAN)
% AS WELL AS THE ROTATION MATRIX SO THAT THE DILATATIONAL PART OF
% THE DEFORMATION IS ZERO (AND IF CANNOT BE RENDER ZERO, THEN THE PROGRAM
% WILL COME TO A HALT, BECAUSE THE CURRENT VERSION (3-March-2023) cannot
% handle accurately distorted domains )
TRANS_COOR = [] ;

% THIS METHOD DOES NOT REQUIRE TO COMPUTE THE JACOBIAN OF THE
% TRANSFORMATION


% X_e = TRANSLATION_ce + DILAT_FACTOR*Q*(Xref) ;
% The preceding expression assumes that Xref are the coordinates of the
% parent domain with respect to its centroid.
% However, such centroid is calculated using the underlying FE mesh. In
% the case of the physical domain, we cannot calculate such centroid
% because there are no underlying FE mesh.
% For instance, if there is neither rotation nor dilation, then
% X_e  -Ce  =    Xref - C_ref ---> X_e = C_e + (Xref-C_ref)
% If we introduce rotations
% X_e = C_e + DILAT_FACTOR*Q*(Xref-C_ref) ;



Xref=  EIFEoper.INFO.FESHAPE_coarse_elem_transf_coord.COOR'; % Coordinates parent domain
CNref = EIFEoper.INFO.FESHAPE_coarse_elem_transf_coord.CN;  % Connectivities parent domain
Xref= Xref(:,CNref) ;
 
C_ref = sum(Xref,2)/size(Xref,2) ; % Pseudo-centroid parent domain
XrefREL = bsxfun(@minus,Xref,C_ref) ; % Relative coordinates parent domain

% Compute maximum length parent domain
% JAHO, 17-APRIL-2024
% It is assumed here that the connectivities are the same for both the
% reference domain and the pyhysical domain
% --------------------------------
METHOD_COMPUTE_maximum_length = 0;  % 0  Before 17-Apr-2024
dLref = DistanceMaximPointsEIFEM(Xref,METHOD_COMPUTE_maximum_length) ;
% -------------------------------------
% Translation. First compuyte the pseudo-centroid of the physical
% domain
C_e = sum(Xe,2)/size(Xe,2) ;
XeREL = bsxfun(@minus,Xe,C_e) ; % Relative coordinates
TRANSLATION_ce = C_e ;  % This is the translation
% --------------------------------------
% Scale factor
%     XAUG = [XeREL,XeREL(:,1)] ;
%     dL  =diff(XAUG,1,2) ;
%     dLE = max(sqrt(sum(dL.^2,1))) ;
dLE = DistanceMaximPointsEIFEM(XeREL,METHOD_COMPUTE_maximum_length) ;
factorSCALE = dLref/dLE ;

% ndim
ndimFINE = size(EIFEoper.MESH.COOR,2) ; 

detJe = 1/(factorSCALE^ndimFINE) ;
TRANS_COOR.ROTATION = 1;
TRANS_COOR.TRANSLATION = TRANSLATION_ce;
TRANS_COOR.CenterRotationREFERENCE = C_ref;
TRANS_COOR.SCALEFACTOR = 1/factorSCALE;
TRANS_COOR.detJe = detJe ;
TRANS_COOR.IndexParentDomain = 1;  





%     else
%         if abs(detJe-detJe_ref)/detJe_ref >= 1e-10
%             error('Only constant Jacobian transformations allowed')
%         end