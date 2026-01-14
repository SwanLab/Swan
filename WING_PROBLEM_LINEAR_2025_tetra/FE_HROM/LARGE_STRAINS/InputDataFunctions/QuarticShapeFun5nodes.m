function  [Nshape,xi,IndPoints] =   QuarticShapeFun5nodes(ElemBnd_loc,COORnodes,COORbnd,Nshape,xi_nodes)
 %--------------------------------------------------------------------------
% function [Nshape, xi, IndPoints] = QuarticShapeFun5nodes(ElemBnd_loc, COORnodes, COORbnd, Nshape, xi_nodes)
%
% PURPOSE:
%   Evaluates Lagrange-type shape functions of arbitrary order on a straight
%   boundary segment defined by a set of nodes. It identifies which boundary
%   nodes lie on this segment, maps them to local coordinates in [-1,1], and 
%   computes the corresponding shape function values using the provided or default
%   reference positions.
%
% INPUTS:
%   - ElemBnd_loc : 1xP vector with indices of P boundary nodes that define the
%                   interpolation segment (ordered from left to right).
%   - COORnodes   : Nnodes x 2 matrix of global coordinates for all mesh nodes.
%   - COORbnd     : Nbndnodes x 2 matrix of coordinates for all boundary nodes.
%   - Nshape      : (Nbndnodes x Ntotalnodes) matrix to be updated with shape
%                   function values for identified boundary points.
%   - xi_nodes    : (Optional) 1xP vector of reference positions (in [-1,1]) for the
%                   interpolation nodes. If empty or not provided, defaults to 
%                   equispaced points in [-1,1].
%
% OUTPUTS:
%   - Nshape      : Updated shape function matrix with contributions for this element.
%   - xi          : Local coordinates in [-1,1] of the identified boundary points.
%   - IndPoints   : Indices in COORbnd of boundary points located on the segment.
%
% METHOD:
%   1. Computes the direction and length of the segment from its first and last node.
%   2. Projects all boundary nodes onto this segment and retains only those that lie
%      close to the segment line and within its bounds (via geometric tolerances).
%   3. Maps valid boundary nodes to the local coordinate system in [-1,1].
%   4. Evaluates Lagrange shape functions at these mapped coordinates.
%   5. Populates the corresponding rows in Nshape.
%
% ASSUMPTIONS:
%   - The element is defined over a straight segment (curved geometry is not supported).
%   - Tolerances are used to decide point inclusion, so very close points are treated as lying on the segment.
%
% NOTES:
%   - This generalizes the original quartic (5-node) implementation to support any number of nodes.
%   - Designed for high-order FEM applications along 1D interfaces embedded in 2D.
%   - Useful for applying boundary or interface conditions in a variational formulation.
%
% USAGE EXAMPLE:
%   xi_nodes = linspace(-1,1,P);  % For P nodes
%   [Nshape_updated, xi_coords, bnd_indices] = QuarticShapeFun5nodes(ElemBnd_loc, COORnodes, COORbnd, Nshape, xi_nodes);
%
% AUTHOR:
%   J.A. Hernández Ortega (JAHO), updated 22-May-2025, Universitat Politècnica de Catalunya (UPC), Terrassa
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
if nargin == 0
    load('tmp1.mat')
end
if isempty(xi_nodes)
  xi_nodes =   [-1, -0.5, 0, 0.5, 1] ; 
end


x1 = COORnodes( ElemBnd_loc(1),:) ;  % First point segment (assuming it is straight)
x2 = COORnodes(ElemBnd_loc(end),:) ;
norm_x1x2 = norm(x2-x1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Compute the coordinates of all the boundary nodes with
% respect to the first point of the interface we are considering
COORbnd_rel = bsxfun(@minus,COORbnd',x1(:))' ;
% Orthogonalize the vectors
norm_COORbnd_rel = sqrt(sum(COORbnd_rel.^2,2)) ;
% See if any of the points is the reference point
TOL_zero  = 1e-6;
[DDD,IIII]=min(norm_COORbnd_rel) ;
CandidatePoints = 1:length(norm_COORbnd_rel) ;
if DDD/norm_x1x2  <TOL_zero
    IndPointsInclude= IIII ;
    CandidatePoints = setdiff(1:length(norm_COORbnd_rel),IndPointsInclude) ;
end
% Normalization
uBND = bsxfun(@times,COORbnd_rel(CandidatePoints,:),1./norm_COORbnd_rel(CandidatePoints)) ;
% Compute scalar product with  the vector that goes from x1 to x2
u  = (x2-x1)/(norm(x2-x1)) ;
PROY = u*uBND' ;
% A necessary condition ---but not sufficient--- is that the nodes should
% lie on the  line defined by x2 and x1
TOL_zero = 1e-4;
IndElements_1 = find(abs(PROY-1)/norm_x1x2<=TOL_zero );
% The other condition is that COORbnd_rel divided by the norm of x2-x1
% should be less than one
Proy_u = u*COORbnd_rel(CandidatePoints(IndElements_1),:)'/norm_x1x2 ;
TOL  = 1e-5 ;
SubIndElements_2 = find(Proy_u <= (1+TOL) & Proy_u>=TOL);
IndPoints = IndElements_1(SubIndElements_2) ;

IndPoints= unique([IndPointsInclude,CandidatePoints(IndPoints)]) ;


xREF = (x1 + x2) / 2; 
COORbnd_rel_0 = 2 * bsxfun(@minus, COORbnd', xREF(:))' / norm_x1x2;
xi = u * COORbnd_rel_0(IndPoints, :)';

% Define reference xi positions for 5 quartic nodes
num_nodes = length(xi_nodes);

% Loop over each node to compute its shape function
for i = 1:num_nodes
    imode = ElemBnd_loc(i);
    N_i = ones(size(xi));
    for j = 1:num_nodes
        if j ~= i
            N_i = N_i .* (xi - xi_nodes(j)) / (xi_nodes(i) - xi_nodes(j));
        end
    end
    Nshape(IndPoints, imode) = N_i;
end