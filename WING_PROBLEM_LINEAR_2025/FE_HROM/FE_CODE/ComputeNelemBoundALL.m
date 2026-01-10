function [ NelemB wSTb ] = ComputeNelemBoundALL(COOR,CNb,TypeElementB) 
%%%%
% -------------------------------------------------------------------------
% COMMENTS (generated automatically by ChatGPT on 7-Nov-2025)
%
% PURPOSE:
%   Compute the shape-function matrix and quadrature weights for all
%   boundary (surface) finite elements in a vectorized manner.
%   The subroutine returns:
%     1) NelemB : (nelemB*ngausB x nnodeEb) matrix containing, stacked by
%                 element and Gauss point, the shape functions (N-matrices)
%                 of all boundary elements.
%     2) wSTb   : (nelemB*ngausB x 1) vector containing the product of Gauss
%                 weights and Jacobian determinants for each boundary Gauss
%                 point across all elements.
%
% INPUTS:
%   COOR          - (nnode x ndim) matrix of nodal coordinates.
%   CNb           - (nelemB x nnodeEb) connectivity matrix for boundary elements.
%   TypeElementB  - string indicating the type of boundary element 
%                   ('Line', 'Triangle', 'Quadrilateral', etc.).
%
% OUTPUTS:
%   NelemB  - Matrix of boundary shape functions for all Gauss points.
%   wSTb    - Vector with the product (detJb * w_g) for all boundary Gauss points.
%
% DETAILS:
%   - The routine assumes that boundary elements live in an (ndim-1)-dimensional
%     parametric space (e.g., 2D boundaries in a 3D domain, or 1D edges in 2D).
%   - Vectorized implementation: all boundary elements are processed simultaneously.
%   - The boundary Jacobian determinant detJb is obtained through:
%        XeALL  = COORallELEM(...)   → elementwise nodal coordinates
%        SeALL  = ChangeCoordBndVect(XeALL, ndimB)  → coordinate projection
%        JeALL  = SeALL * (dN/dξ)ᵀ   → Jacobian matrices (ndimB x ndimB)
%        detJeALL = determinantVECTORIZE(JeALL, ndimB)
%   - The Gauss quadrature weights and shape functions are obtained from:
%        [weig, posgp, shapef, dershapef] =
%             ComputeElementShapeFun(TypeElementB, nnodeEb, 'RHS')
%   - The matrix of shape functions NelemB is constructed by tiling shapef
%     for all boundary elements, since shape functions are identical for all
%     elements in the reference space.
%
% COMPUTATION OF wSTb:
%   - For each Gauss point g, detJeALL (size nelemB x 1) is multiplied by the
%     corresponding Gauss weight w_g, replicated for each element:
%         wLOCa = detJeALL .* weigREP(g:ngaus:nelemB*ngaus)
%     and stored in the positions corresponding to Gauss point g across all
%     boundary elements:
%         wSTb(g:ngaus:nelemB*ngaus) = wLOCa
%
% SANITY / VALIDATION:
%   - If ChangeCoordBndVect returns NaNs, an error is raised to indicate an
%     invalid local coordinate transformation.
%   - Negative Jacobians are not filtered here but could indicate inverted
%     boundary element orientation.
%
% AUTHOR / HISTORY:
%   Joaquín A. Hernández (jhortega@cimne.upc.edu), 27-Oct-2015
%   Comments clarification (this header): 7-Nov-2025
% -------------------------------------------------------------------------



%%%%
% This subroutine   returns
%  1) the matrix Nelem (ndim*nelemB*ngausB x nnodeE)  consisting of all boundary element N-matrices
%  2) wSTb = zeros(nelemB*ngausB,1) ;  % Vector containinig the product of weights and Jacobians at all boundary gauss points
% Inputs:   COOR: Coordinate matrix (nnode x ndim), % CNb: Connectivity
% matrix  for boundary elements
% TypeElement: Type of finite element (quadrilateral,...),
%% Vectorized version
% Joaquín A. Hernández (jhortega@cimne.upc.edu), 27-Oct-2015
%dbstop('12')
if nargin == 0
    load('tmp1.mat')
end
nnode = size(COOR,1); ndim = size(COOR,2)  ; ndimB= ndim -1; nelemB = size(CNb,1); nnodeEb = size(CNb,2) ;
% nstrain = size(celasglo,1) ;
% Shape function routines (for calculating shape functions and derivatives)
TypeIntegrand = 'RHS';
[weig,posgp,shapef,dershapef] = ComputeElementShapeFun(TypeElementB,nnodeEb,TypeIntegrand) ;
ngaus = length(weig) ;
% Therefore, the sought-after Nelem matrix can be computed by simply making  nelem tiling
% copies of Ne
NelemB = repmat(shapef,nelemB,1);
% -------------------------------
%%%% Computation of the boundary weights wSTb
% -------------------------------
wSTb = zeros(nelemB*ngaus,1) ;  % Vector containinig the product of weights and Jacobians at all gauss points
% COORDINATE MATRIX (for boundary elements) arranged in a nelemB*ndim x nnodeEB matrix
XeALL= COORallELEM(ndim,nelemB,nnodeEb,CNb,COOR) ;
% Change of coordinates
SeALL = ChangeCoordBndVect(XeALL,ndimB) ;
if any(isnan(SeALL))
    error('Some issue has poped out when trying to make the coordinates change')
end
% Let us define a matrix ROWSgauss such that   Belem(ROWSgauss(g,:),:) returns
% the B-matrices of the g-th points of all elements
weigREP = repmat(weig',nelemB,1)  ; % nelem x 1 tiling copies of weig

for  g = 1:ngaus
    % Matrix of derivatives for Gauss point "g"
    BeXi = dershapef(:,:,g) ;
    % ----------------------------------------
    % Coordinate change (vectorized form)
    
    
    %%%%%
    % Jacobian Matrix for the g-th G. point of all elements %
    JeALL = SeALL*BeXi' ;
    %%%%%%%%%
    % JAcobian
    detJeALL= determinantVECTORIZE(JeALL,ndimB) ;
    % Weight vectors
    % --------------
    wLOCa = detJeALL.*weigREP(g:ngaus:nelemB*ngaus) ;
    wSTb(g:ngaus:nelemB*ngaus) =  wLOCa ;
end

end
