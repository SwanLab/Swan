function [BNDINTERFACE,CNbLOC_BNDINFERFACE_ALL] = BoundaryINFOassign(InterfaceBoundariesLocalNumbering,MESH,OPERFE)
%--------------------------------------------------------------------------
% FUNCTION: BoundaryINFOassign
%
% PURPOSE:
%   Retrieves geometric and finite element information for the specified set
%   of interface boundaries (faces) in a given subdomain mesh. It builds a
%   structured representation (`BNDINTERFACE`) of each interface boundary
%   including connectivity, nodes, centroid, shape functions, normals, and
%   Gauss integration data, which are required for the definition of 
%   interface modes or assembly of interface-related operators.
%
% INPUTS:
%   - InterfaceBoundariesLocalNumbering : Array of indices specifying the local
%     numbering of the boundaries in MESH that are considered "interfaces".
%
%   - MESH : Structure containing mesh data, including:
%       * .Indexes_faces_bnd_element → Cell array with element indices per face
%       * .CNb                       → Boundary connectivity matrices per face
%       * .PROPERTIES_FACES         → Face properties (area, centroid, shape functions, etc.)
%
%   - OPERFE : Structure with precomputed finite element operators, including:
%       * .NstT_W_N_boundaries      → Assembled left-weighted N-matrices per face
%
% OUTPUTS:
%   - BNDINTERFACE : Array of structs (1×nFaces) with fields:
%       * .CNb                         → Local connectivities of each boundary
%       * .NODES                       → Unique nodes of the boundary
%       * .AREA                        → Area (or length in 2D) of the boundary
%       * .COORrelA_global            → Coordinates relative to boundary centroid
%       * .CENTROID                   → Global centroid coordinates of the face
%       * .GeometricMassMatrix        → Mass matrix computed over the boundary
%       * .Nst                        → Shape function matrix at Gauss points
%       * .wST                        → Gauss weights over the boundary
%       * .UnitNormalAtGaussPoint     → Normal vectors at Gauss points
%       * .UnitTangentAtGaussPoint    → Tangent vectors at Gauss points
%       * .NstT_W_N_boundaries        → Mass-like term: N' * W * N
%
%   - CNbLOC_BNDINFERFACE_ALL : Global connectivity matrix formed by vertically
%     stacking all interface boundary connectivities.
%
% APPLICATION:
%   This function is central in domain decomposition, reduced-order modeling
%   (e.g., EIFEM), and multi-scale methods, where interface representation
%   and interaction with neighboring subdomains are needed.
%
% EXAMPLE USAGE:
%   [BNDinterface,CNbLOC] = BoundaryINFOassign([1,3],MESH,OPERFE);
%
% SEE ALSO:
%   GeometricVarDOMAINS, GeometricVarDOMAINScINTF, ConstructBasisRigidBody_MASSM
%
% AUTHOR:
%   Joaquín A. Hernández Ortega (CIMNE/UPC)
%   Last modified: 28-Jan-2025
%   Comments by ChatGPT4, 13-May-2025
%--------------------------------------------------------------------------



% see /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/DOMAINdecompHROM/MultiHROM/GeometricVarDOMAINS.m

BNDINTERFACE =[] ;% cell(length(InterfaceBoundariesLocalNumbering),1) ;

CNbLOC_BNDINFERFACE_ALL = cell(length(InterfaceBoundariesLocalNumbering),1) ; 

for iface = 1:length(InterfaceBoundariesLocalNumbering)
    ilocFACE = InterfaceBoundariesLocalNumbering(iface) ;
    IndexElemBndLoc = MESH.Indexes_faces_bnd_element{ilocFACE} ;
    CNbLOC = MESH.CNb(IndexElemBndLoc,:)  ;
    CNbLOC_BNDINFERFACE_ALL{iface} = CNbLOC ; 
    NODES = unique(CNbLOC(:)) ;
    BNDINTERFACE(iface).CNb = CNbLOC;  % Connectivities
    BNDINTERFACE(iface).NODES = NODES;  % Nodes, sorted in ascending order
    BNDINTERFACE(iface).AREA = MESH.PROPERTIES_FACES{ilocFACE}.AREA ;  % Area/length
    BNDINTERFACE(iface).COORrelA_global = MESH.PROPERTIES_FACES{ilocFACE}.COORrelA_global ;  % Coordinates relative to its centroid
    BNDINTERFACE(iface).CENTROID = MESH.PROPERTIES_FACES{ilocFACE}.CENTROID ; % Centroid (referred to the centroid of the domain)
    BNDINTERFACE(iface).GeometricMassMatrix = MESH.PROPERTIES_FACES{ilocFACE}.GeometricMassMatrix ;  % Geometric mass matrix
    BNDINTERFACE(iface).Nst = MESH.PROPERTIES_FACES{ilocFACE}.Nst ;  % displacement_at_gauss_points = Nst*displacement_at_nodes
    BNDINTERFACE(iface).wST = MESH.PROPERTIES_FACES{ilocFACE}.wST ;  % Gauss weights
    BNDINTERFACE(iface).UnitNormalAtGaussPoint = MESH.PROPERTIES_FACES{ilocFACE}.UnitNormalAtGaussPoint ;  % Normals
    BNDINTERFACE(iface).UnitTangentAtGaussPoint = MESH.PROPERTIES_FACES{ilocFACE}.UnitTangentAtGaussPoint ;  % Normals
     BNDINTERFACE(iface).NstT_W_N_boundaries = OPERFE.NstT_W_N_boundaries{ilocFACE}  ;  % Normals
end

CNbLOC_BNDINFERFACE_ALL = cell2mat(CNbLOC_BNDINFERFACE_ALL) ; % 28-Jan-2025