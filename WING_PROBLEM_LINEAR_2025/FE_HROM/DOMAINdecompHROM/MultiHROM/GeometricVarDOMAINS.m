function [faceDOFS,PhiRB,PsiRBf,Mdom,MdomFF,MESH] =   GeometricVarDOMAINS(OPERFE,MESH,DATA,DATAcommon)
%--------------------------------------------------------------------------
% FUNCTION: GeometricVarDOMAINS
%
% PURPOSE:
%   Computes geometric and kinematic quantities associated with a subdomain
%   used in multiscale or reduced-order finite element methods (e.g., EIFEM).
%   The function identifies the interface boundaries of the domain, extracts 
%   associated degrees of freedom (DOFs), and computes projection operators 
%   based on rigid-body kinematics.
%
%   The local coordinate system is assumed to be centered at the **centroid 
%   of the domain** (not the interface). For interface-centered coordinate
%   definitions, use `GeometricVarDOMAINScINTF`.
%
% USAGE:
%   [faceDOFS, PhiRB, PsiRBf, Mdom, MdomFF, MESH] = ...
%       GeometricVarDOMAINS(OPERFE, MESH, DATA, DATAcommon)
%
% INPUT:
%   - OPERFE     : Structure containing finite element operators (B-matrix, weights, etc.)
%   - MESH       : Mesh structure with boundary faces, node coordinates, etc.
%   - DATA       : Structure with general model information (e.g., dimension, number of DOFs)
%   - DATAcommon : Auxiliary configuration options, including:
%                  * InterfaceBoundariesLocalNumbering : IDs of the interface boundaries
%                  * LateralBoundariesLocalNumbering   : IDs of lateral (non-coupled) boundaries
%
% OUTPUT:
%   - faceDOFS   : Cell array containing DOFs on each interface face
%   - PhiRB      : Rigid-body modes of the full subdomain (w.r.t. centroid)
%   - PsiRBf     : Projection of rigid-body modes onto interface boundaries via mass matrix
%   - Mdom       : Geometric (scalar) mass matrix of the domain (used for modal projections)
%   - MdomFF     : Block-diagonal geometric mass matrix for the interface (vector-valued)
%   - MESH       : Updated mesh structure with added fields:
%                   * BNDinterface, BNDlateral (boundary face definitions)
%                   * faceDOFSall, faceNODESall (aggregated interface boundary info)
%
% KEY STEPS:
%   1. Identify interface and lateral boundaries based on provided labels.
%   2. Extract DOFs and nodes associated with the interface boundaries.
%   3. Recover domain-level rigid-body modes (`PhiRB`) precomputed with respect to centroid.
%   4. Compute geometric mass matrix of the domain (`Mdom`).
%   5. Assemble geometric interface mass matrix (`MdomFF`) in vector format.
%   6. Project rigid-body modes onto the interface using `PsiRBf = MdomFF * PhiRB`.
%
% APPLICATIONS:
%   - Projection of boundary kinematics in empirical interscale methods (EIFEM)
%   - Interface reduction in component-based ROMs
%   - Modal matching and enforcement of rigid compatibility across interfaces
%
% REMARKS:
%   - Interface rigid-body coupling is evaluated using `Grb = PhiRB^T * PsiRBf`
%   - If needed, the function `GeometricVarDOMAINScINTF` can be used to align
%     the coordinate system with the interface boundary centroid.
%
% SEE ALSO:
%   - GeometricVarDOMAINScINTF
%   - BoundaryINFOassign
%   - CompWeightDiag
%   - RigidBodyModes_fromCOOR
%
% AUTHOR:
%   Joaquín A. Hernández Ortega
%   UPC - CIMNE, Barcelona
%   Date: 28-Jan-2025
%   Comments by ChatGPT4, 13-May-2025
%--------------------------------------------------------------------------


% This function returns, for the studied subdomain whose mesh and finite
% element operators are given in OPERFE and MESH
%  faceDOFS:

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
MESH.BNDinterface = BoundaryINFOassign(DATAcommon.InterfaceBoundariesLocalNumbering,MESH,OPERFE) ;




% -----------------------------------------------------------

% Lateral interfaces. Retrieving information
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

 

% --------------------------------------------------

% ------------------------------------------------
% Rigid-body modes ------------------------------
% ----------------------------------------------------
% Computed in /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/FE_CODE/ConstructBasisRigidBody_MASSM.m
% (With respect to the centroid of the domain )
PhiRB  = MESH.GEOproperties.RIGID_BODY_MODES ; 

% -------------------------------------------------------
% Geometric mass matrix 
% -------------------------------------------------------
W =CompWeightDiag(OPERFE.wSTs_RHS,DATA.MESH.ndim) ; 
Mdom = OPERFE.Nst'*(W*OPERFE.Nst) ; 

% Boundary DOFs
% ------------------------
 faceDOFS = cell(1,length(MESH.BNDinterface)) ; 
  faceNODES = cell(1,length(MESH.BNDinterface)) ; 

for iface =1:length(MESH.BNDinterface)
    nodesloc = MESH.BNDinterface(iface).NODES ; 
    MESH.BNDinterface(iface).DOFs = small2large(MESH.BNDinterface(iface).NODES,DATA.MESH.ndim) ; 
    faceDOFS{iface} = MESH.BNDinterface(iface).DOFs ; 
    faceNODES{iface} = nodesloc ; 
end


%----------------------------------------------------------
% Boundary matrix MdomFF
% Boundary nodes are assumed to be ordered in ascending order 
% --------------------------------------------------------
f = unique(cell2mat(faceNODES(:))) ;
ndimLOC = 1;
wDIAGb = CompWeightDiag(OPERFE.wSTb,ndimLOC)  ;
Mbound = OPERFE.NstB_left'*(wDIAGb*OPERFE.NstB_left) ;
Mff = Mbound(f,f) ;
% ------------------------------------------------------------- 
% This is a matrix relating scalar nodal vectors. Extension to vectors of
% ndim entries per node
nnodeB = size(Mff,1) ;
ndim = DATA.MESH.ndim ;
MdomFF = sparse(nnodeB*DATA.MESH.ndim,nnodeB*ndim) ;
for idim = 1:DATA.MESH.ndim
    MdomFF(idim:ndim:end,idim:ndim:end)  = Mff ;
end
faceDOFSall = small2large(f,DATA.MESH.ndim) ; 
PsiRBf = MdomFF*PhiRB(faceDOFSall,:) ; 

MESH.faceDOFSall = faceDOFSall ; 
MESH.faceNODESall = f ;  