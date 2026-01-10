function [CentroidFA,AREA,Mst,Nst,wSTb] =CentroidGeometricMassMatrixNEW(COOR,nodesfA,CONNECTb,TypeElementB)
%--------------------------------------------------------------------------
% FUNCTION: CentroidGeometricMassMatrixNEW
%
% PURPOSE:
%   Computes the geometric mass matrix, shape function matrix, and centroid
%   of a boundary face composed of one or more boundary elements.
%
%   This function is used primarily in the context of multiscale methods
%   (e.g., EIFEM) where mass lumping or projection operations are required
%   on the boundary of a subdomain (interface coupling).
%
% USAGE:
%   [CentroidFA, AREA, Mst, Nst, wSTb] = ...
%       CentroidGeometricMassMatrixNEW(COOR, nodesfA, CONNECTb, TypeElementB)
%
% INPUT:
%   - COOR         : Nodal coordinates of the full domain (size: nnode x ndim)
%   - nodesfA      : List of nodes on the current boundary face (can be empty; if so, inferred from CONNECTb)
%   - CONNECTb     : Connectivity of boundary elements (rows = elements, cols = nodes per element)
%   - TypeElementB : Type of boundary element (e.g., 'LIN2', 'QUAD4', etc.)
%
% OUTPUT:
%   - CentroidFA   : Coordinates of the geometric centroid of the face
%   - AREA         : Total area (or length in 2D) of the boundary face
%   - Mst          : Geometric mass matrix (lumped or consistent)
%   - Nst          : Global shape function matrix (evaluated at Gauss points)
%   - wSTb         : Vector of quadrature weights associated with boundary Gauss points
%
% KEY STEPS:
%   1. If `nodesfA` is empty, it is reconstructed from `CONNECTb`.
%   2. The connectivity `CONNECTb` is renumbered locally.
%   3. Gauss points and weights are computed using `ComputeNelemBoundALL`.
%   4. The global shape function matrix `Nst` is assembled for the boundary.
%   5. A diagonal weight matrix `wSTdiag` is built and used to compute the
%      consistent mass matrix `Mst = Nstᵀ W Nst`.
%   6. The geometric centroid is computed by weighted averaging of nodal coordinates.
%
% NOTES:
%   - This function assumes scalar fields (ndim = 1) for mass integration.
%   - `Mst` is often used to project quantities onto the boundary or define interface metrics.
%   - If quadrature weights contain `NaN`, a connectivity renumbering error is triggered.
%
% SEE ALSO:
%   - ComputeNelemBoundALL
%   - AssemblyNGlobal
%   - CompWeightDiag
%   - RenumberConnectivities
%
% AUTHOR:
%   Joaquín A. Hernández Ortega
%   UPC - CIMNE, Barcelona
%   Date: 22-Feb-2025
%   Comments by ChatGPT4, 13-May-2025
%--------------------------------------------------------------------------

if nargin == 0
    load('tmp3.mat')
end

if isempty(nodesfA)
    nodesfA = unique(CONNECTb(:)) ; 
end
nodesfA = sort(nodesfA) ; 
COOR_FACE = COOR(nodesfA,:) ; % Coordinates of this face

CNface1= RenumberConnectivities(CONNECTb  ,1:length(nodesfA)) ;


[ Nelem, wSTb ] = ComputeNelemBoundALL(COOR_FACE,CNface1,TypeElementB) ; 
if any(isnan(wSTb))
    error('Error in renumbering connectivities')
end
% -------------------------------------------------------------------------
% Assembly 
nelem = size(CNface1,1) ; nnodeE = size(CNface1,2) ; ndim = 1; 
ngaus = size(Nelem,1)/nelem ; nnode = size(nodesfA,1) ;  
Nst = AssemblyNGlobal(Nelem,nelem,nnodeE,ndim,ngaus,CNface1,nnode) ;

 
wSTdiag = CompWeightDiag(wSTb,1)  ;
Mst = (wSTdiag*Nst)'*Nst ;
% Recomputing centroid
CentroidFA = zeros(1,size(COOR_FACE,2)) ;
AREA = sum(wSTb) ;
for idim = 1:size(COOR_FACE,2)
    CentroidFA(idim) = wSTb'*(Nst*COOR_FACE(:,idim))/AREA ;
end