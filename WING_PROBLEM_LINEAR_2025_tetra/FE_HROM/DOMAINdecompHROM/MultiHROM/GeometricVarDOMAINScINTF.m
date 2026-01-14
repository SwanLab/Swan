function [PhiRB,PsiRBf,Mdom,MdomFF,MESH] =   GeometricVarDOMAINScINTF(OPERFE,MESH,DATA,DATAcommon)
 %--------------------------------------------------------------------------
% FUNCTION: GeometricVarDOMAINScINTF
%
% PURPOSE:
%   Computes geometric and kinematic variables for a finite element subdomain 
%   in the context of multiscale or reduced-order modeling. Unlike its 
%   companion function `GeometricVarDOMAINS`, this version defines the 
%   local coordinate system origin at the **centroid of the interface boundary**, 
%   not the domain centroid. This is critical for problems where interface 
%   behavior dominates, such as cohesive contact, beam-like elements, or 
%   embedded interfaces in ROMs.
%
% USAGE:
%   [PhiRB, PsiRBf, Mdom, MdomFF, MESH] = ...
%       GeometricVarDOMAINScINTF(OPERFE, MESH, DATA, DATAcommon)
%
% INPUT:
%   - OPERFE     : Structure containing FE operators (e.g., shape functions, weights).
%   - MESH       : Mesh structure with domain and boundary connectivity.
%   - DATA       : General problem data (e.g., number of DOFs, dimensions).
%   - DATAcommon : Auxiliary input options, especially:
%                   * InterfaceBoundariesLocalNumbering (indexes of interface faces)
%                   * LateralBoundariesLocalNumbering (optional)
%
% OUTPUT:
%   - PhiRB      : Rigid-body modes for the subdomain, defined with respect to the 
%                  interface centroid.
%   - PsiRBf     : Projection of rigid-body modes onto the interface boundary, 
%                  via geometric mass matrix.
%   - Mdom       : Geometric (consistent) mass matrix of the full domain.
%   - MdomFF     : Geometric mass matrix of the interface (block-diagonal by component).
%   - MESH       : Updated mesh structure, with fields like:
%                   * COOR (shifted w.r.t interface centroid),
%                   * faceNODESall, faceDOFSall (interface identification),
%                   * BNDinterface, BNDlateral (boundary topology),
%                   * GEOproperties.CENTROID (updated),
%                   * GEOproperties.RIGID_BODY_MODES (w.r.t interface)
%
% KEY STEPS:
%   1. Extracts connectivity and node information for interface boundaries.
%   2. Computes centroid and geometric mass matrix of the interface (CentroidGeometricMassMatrixNEW).
%   3. Shifts the mesh coordinates so the interface centroid becomes the origin.
%   4. Recomputes rigid-body modes (`PhiRB`) and their projection onto interface DOFs (`PsiRBf`).
%   5. Computes the domain-level geometric mass matrix (`Mdom`) and interface version (`MdomFF`).
%   6. Returns updated mesh and interface-local kinematics.
%
% APPLICATIONS:
%   - Corotational ROMs (e.g., EIFE/ECM in large rotations)
%   - Substructuring with non-matching interfaces
%   - Detection of coupling terms between rigid-body translations and rotations
%     via inspection of matrix `Grb = PhiRB^T * PsiRBf`
%
% SEE ALSO:
%   - GeometricVarDOMAINS
%   - CentroidGeometricMassMatrixNEW
%   - RigidBodyModes_fromCOOR
%   - BoundaryINFOassign
%
% AUTHOR:
%   Joaquín A. Hernández Ortega
%   UPC-CIMNE, Barcelona
%   Date: 28-Jan-2025
%   Comments by ChatGPT4, 13-May-2025
%--------------------------------------------------------------------------
 
 




% This function returns, for the studied subdomain whose mesh and finite
% element operators are given in OPERFE and MESH
%  Similar to GeometricVarDOMAINS.m, but with the origin of the local axes
%  on the centroid of the boundary interface, rather than the centroid of
%  the domain 
% JAHO, 28-JAN-2025, Balmes 185, Barcelona 
% -------------------------------------------------------------------------

if nargin == 0
    load('tmp1.mat')
end

% 1) Number of boundary interfaces
% -----------------------------------------
% Total number of boundary surfaces
NumberBoundSurfaces = length(MESH.PROPERTIES_FACES) ;

% For instannce, for a 2D quadrilateral subdomain:
% MESH.PROPERTIES_FACES
%
% ans =
%
%   1×4 cell array
%
%     {1×1 struct}    {1×1 struct}    {1×1 struct}    {1×1 struct}

% OUT OF THIS 4 SURFACES, WE HAVE TO SPECIFY WHICH ARE LATERAL BOUNDARIES,
% AND INTERFACE BOUNDARIES

DATAcommon = DefaultField(DATAcommon,'InterfaceBoundariesLocalNumbering',1:NumberBoundSurfaces) ;

% By default, all boundaries are regarded as interface boundaries. In this
% case, for instance, if we want to turn this into a "beam" element, we
% should select only 2 of this boundaries, say
% DATA.InterfaceBoundariesLocalNumbering  = [1,3]


% --------------------
% Boundary interfaces. Retrieving information 
%---------------------
[MESH.BNDinterface,CN_BNDinterface_all] = BoundaryINFOassign(DATAcommon.InterfaceBoundariesLocalNumbering,MESH,OPERFE) ;

% Computing centroid, area, and mass matrix (scalar) INTERFACE BOUNDARY AS
% A WHOLE 
 [Centroid_interface,AREA_interface,MffNEW,Nst_intf,wSTb_intf] =CentroidGeometricMassMatrixNEW(MESH.COOR,[],CN_BNDinterface_all,MESH.TypeElementB) ;
% RECOMPUTING COORDINATES MESH.COOR 
disp('Recomputing  coordinates domain (with respect to centroid interface boundary)')
for idim = 1:length(Centroid_interface)
    MESH.COOR(:,idim) = MESH.COOR(:,idim) - Centroid_interface(idim); 
end
MESH.COOR_Centroid_interface_respect_CENTROID = Centroid_interface ; 
% RECOMPUTING RIGID BODY MODES with respect centrod interfaces 
PhiRB = RigidBodyModes_fromCOOR(MESH.COOR) ; 

% REDEFINING RIGID BODY MODES AND CENTROID 
disp('Redefining the rigid boyd modes/location of the origin of the local coordinate system (= centroid interface boundary)')
 MESH.GEOproperties.RIGID_BODY_MODES = PhiRB ; 
 MESH.GEOproperties.CENTROID = MESH.GEOproperties.CENTROID +  Centroid_interface ; 


% lIST OF BOUNDARY NODES interface boundaries
f = unique(CN_BNDinterface_all) ; 
 

% INTERFACE GEOMETRIC MASS MATRIX, FOR ALL SPATIAL DIMENSIONS
nnodeB = size(MffNEW,1) ;

ndim = DATA.MESH.ndim ; %
MdomFF = sparse(nnodeB*DATA.MESH.ndim,nnodeB*ndim) ;
for idim = 1:DATA.MESH.ndim
    MdomFF(idim:ndim:end,idim:ndim:end)  = MffNEW ;
end
faceDOFSall = small2large(f,DATA.MESH.ndim) ;  % DOFs associated to INTERFACE BOUNDARY NODES 



PsiRBf = MdomFF*PhiRB(faceDOFSall,:) ; 

MESH.faceDOFSall = faceDOFSall ; 
MESH.faceNODESall = f ; 

Grb=  PhiRB(faceDOFSall,:)'*PsiRBf ; 

% -----------------------------------
if  DATA.MESH.ndim  == 2 
    GrbTRrt = Grb(1:2,3); 
    
else
    GrbTRrt = Grb(1:3,4:6) ; 
end
Gred  = norm(GrbTRrt)/Grb(1,1); 
disp('Since coordinates are referred to the centroid of the interface')
disp(['norm(GrbTRrt)/Grb(1,1) = ',num2str(Gred)])
disp(['(Translations/rotations interf. boundary are uncoupled)'])

 
 
 
 
% -----------------------------------------------------------

% Lateral interfaces (non-interface boundaries). Retrieving information
% -------------------------------------------------
LateralBoundariesLocalNumbering = setdiff(1:NumberBoundSurfaces,DATAcommon.InterfaceBoundariesLocalNumbering) ;
DATAcommon = DefaultField(DATAcommon,'LateralBoundariesLocalNumbering',LateralBoundariesLocalNumbering) ;
MESH.BNDlateral  = []; 
if ~isempty(DATAcommon.LateralBoundariesLocalNumbering)
    MESH.BNDlateral = BoundaryINFOassign(DATAcommon.LateralBoundariesLocalNumbering,MESH,OPERFE) ;
end
 
MESH = rmfield(MESH,'Indexes_faces_bnd_element') ; 
MESH = rmfield(MESH,'NODES_FACES') ; 
MESH = rmfield(MESH,'PROPERTIES_FACES') ; 


% 
% % --------------------------------------------------
% 
% % ------------------------------------------------
% % Rigid-body modes ------------------------------
% % ----------------------------------------------------
% % Computed in /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/FE_CODE/ConstructBasisRigidBody_MASSM.m
% % (With respect to the centroid of the domain )
% PhiRB  = MESH.GEOproperties.RIGID_BODY_MODES ; 

% -------------------------------------------------------
% Geometric mass matrix 
% -------------------------------------------------------
W =CompWeightDiag(OPERFE.wSTs_RHS,DATA.MESH.ndim) ; 
Mdom = OPERFE.Nst'*(W*OPERFE.Nst) ; 

% 
% faceDOFSall = small2large(f,DATA.MESH.ndim) ; 
% PsiRBf = MdomFF*PhiRB(faceDOFSall,:) ; 
% 





 