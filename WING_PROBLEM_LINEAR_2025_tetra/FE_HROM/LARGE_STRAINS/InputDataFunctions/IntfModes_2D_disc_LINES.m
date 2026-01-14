function [Nshape,COORnodes,OTHER_OUTPUT] =  IntfModes_2D_disc_LINES(INFO_INTERFACE_MESH,CENTROID,COORbnd)
%--------------------------------------------------------------------------
% function [Nshape, COORnodes, OTHER_OUTPUT] = IntfModes_2D_disc_LINES(INFO_INTERFACE_MESH, CENTROID, COORbnd)
%
% PURPOSE:
%   Constructs the shape function matrix `Nshape` for interface modes 
%   defined on a 2D discretized line-based boundary (i.e., interface) as 
%   used in EIFEM (Empirical Interscale Finite Element Method). 
%   The shape functions are computed for mixed interpolation (e.g., linear, 
%   quartic), assuming a 2D geometry, and the interface consists of one or 
%   more line segments defined in `ElemBnd`.
%
%   This function:
%     - Re-centers node coordinates with respect to a given CENTROID.
%     - Computes local shape function values using predefined interpolation 
%       (via ShapeFun_n_nodes) over each line segment.
%     - Assembles these local contributions into a global shape matrix.
%     - Applies a transformation matrix to ensure continuity and proper 
%       weighting of shared interface nodes.
%
% INPUTS:
%   - INFO_INTERFACE_MESH : Structure containing:
%       - COOR : Coordinates of interface nodes (nNodes x 2).
%       - LINES : Cell array defining the node connectivity of each 
%                 boundary segment.
%       - xi_nodes : (Optional) Local reference coordinates on each segment.
%   - CENTROID : 1 x 2 vector specifying the reference point used to 
%                re-center node coordinates.
%   - COORbnd : Coordinates of the boundary points where shape functions 
%               should be evaluated (nPoints x 2).
%
% OUTPUTS:
%   - Nshape : Matrix of shape function values at COORbnd locations 
%                   (nPoints x nModes). It accounts for repeated nodes at 
%                   interface boundaries.
%   - COORnodes : Re-centered coordinates of interface nodes.
%   - OTHER_OUTPUT : Struct with the following fields:
%       - xiEDGES : Local reference coordinates per boundary segment.
%       - IndPointsBNDedge : Indices of boundary points per edge.
%       - NumberNodesPerEdge : Number of nodes per edge (segment).
%
% USAGE NOTES:
%   - This function supports arbitrary numbers of nodes per segment.
%   - Shared nodes between consecutive line segments are weighted by 0.5 
%     (assuming consistent orientation and node ordering).
%   - See reference test file for usage: 
%     /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/...
%     TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/12_REPETITIVE_TRAIN.mlx
%
% AUTHOR:
%   Joaquín A. Hernández Ortega, 26-MAY-2025, Perro de Pavlov, Madrid.  
%   Commented by ChatGPT5
%--------------------------------------------------------------------------

% --------------------------------------------------------------------------------------
if nargin == 0
    load('tmp1.mat')
end
OTHER_OUTPUT = [] ;

INFO_INTERFACE_MESH = DefaultField(INFO_INTERFACE_MESH,'xi_nodes',[]) ;
xi_nodes = INFO_INTERFACE_MESH.xi_nodes;



COORnodes = INFO_INTERFACE_MESH.COOR; %(:,2:end) ;
ElemBnd= INFO_INTERFACE_MESH.LINES ;
nNODES_with_repetition = sum(cellfun(@numel, ElemBnd));

nmodes = size(COORnodes,1) ; % Number of modes
npoints = size(COORbnd,1) ;



Nshape = zeros(npoints,nmodes) ;

% cOORDINATES REFERRED TO THE CENTROID
ndim = size(CENTROID,2);
for idim = 1:ndim
    COORnodes(:,idim) = COORnodes(:,idim) - CENTROID(idim) ;
end

xiEDGES = cell(size(ElemBnd)) ;
IndPointsBNDedge = cell(size(ElemBnd)) ;
NumberNodesPerEdge = zeros(size(ElemBnd)) ;


for ielem = 1:length(ElemBnd)
    
    
    
    xi_nodes = linspace(-1,+1,length(ElemBnd{ielem})) ;
    
    [Nshape,xiEDGES{ielem},IndPointsBNDedge{ielem}] =   ShapeFun_n_nodes(ElemBnd{ielem},COORnodes,COORbnd,Nshape,xi_nodes) ;
    
    NumberNodesPerEdge(ielem) = length(ElemBnd{ielem});
    
    
    
end
OTHER_OUTPUT.xiEDGES = xiEDGES ;
OTHER_OUTPUT.IndPointsBNDedge = IndPointsBNDedge ;
OTHER_OUTPUT.NumberNodesPerEdge = NumberNodesPerEdge ;
 