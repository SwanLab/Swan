function [COORquad, LINES] = generateInterpEdges_GEN2D(COOR,CN, DEGREE_Poly)
%--------------------------------------------------------------------------
% function [COORquad, LINES] = generateInterpEdges_GEN2D(COOR, CN, DEGREE_Poly)
%
% PURPOSE:
%   Generates a set of interpolation points along a list of 2-node edges
%   (segments) in 2D, with the number of interpolation points per edge
%   determined by the corresponding polynomial degree. This is useful, for
%   example, for constructing higher-order finite element geometries or
%   interpolating quantities along interfaces in generalized FEM approaches.
%
% INPUTS:
%   - COOR         : Nnodes x 2 array with coordinates of all original nodes.
%   - CN           : Nelements x 2 array, each row contains the node indices
%                    defining one edge segment.
%   - DEGREE_Poly  : Scalar or Nelements x 1 array with the degree of the
%                    polynomial interpolation for each edge. Each edge will
%                    have (degree + 1) interpolation points. If scalar, it
%                    is applied uniformly to all edges.
%
% OUTPUTS:
%   - COORquad     : Ninterp x 2 array of coordinates for all interpolated
%                    points along all edges. Interpolated points may reuse
%                    endpoints where edges are contiguous.
%   - LINES        : Cell array of size Nelements x 1. Each entry is a vector
%                    of indices into COORquad that corresponds to the local
%                    interpolation points for the corresponding edge.
%
% REMARKS:
%   - The function ensures that interpolation nodes are shared at endpoints
%     for contiguous edges (e.g., when CN(i-1,2) == CN(i,1)).
%   - The last segment is treated cyclically, closing the loop by reusing
%     the initial node index.
%   - Intended for use in 2D mesh preprocessing or interface discretization
%     for generalized/extended FEM methods.
%
% AUTHOR:
%   J.A. Hernández Ortega (UPC), 27-May-2025, Balmes 185, Barcelona
% Comments by ChatGPT4
%--------------------------------------------------------------------------


% Initialize
if nargin == 0
    load('tmp1.mat')
end
COORquad = [];
LINES = cell(size(CN,1), 1);
nodeCount = 0;

nNODES_ALL = 0 ;

if length(DEGREE_Poly) ==1
    DEGREE_Poly = DEGREE_Poly*ones(size(CN,1)) ;
end

for i = 1:length(LINES)
    p = DEGREE_Poly(i) ;
    nPoints = p + 1;
    
    idx1 = CN(i,1);  % Start point index
    idx2 = CN(i,2);      % End point index
    
    P1 = COOR(idx1, :);
    P2 = COOR(idx2, :);
    
    
    
    t = linspace(0, 1, nPoints)';
    edgePoints = (1 - t) * P1 + t * P2;
    
    if   i==1  || (i>1 && (i ~= size(CN,1)) && (CN(i-1,2) ~= CN(i,1)))
        % Either first segment, or a subsequent segment, different from the
        % last one whose first node does not coincide with the last node
        % of the previos one
        COORquad = [COORquad; edgePoints];
        LINES{i} = nodeCount + (1:nPoints);
        nodeCount = nodeCount + nPoints;
    elseif i>1 && (i ~= size(CN,1)) && CN(i-1,2) == CN(i,1)
        % Second or subsequent segments, different from the last one,
        % whose first node coincides with the last one of the previous
        % segment
        COORquad = [COORquad; edgePoints(2:end,:)];
        LINES{i} = (nodeCount-1) + (1:nPoints);
        nodeCount = nodeCount + nPoints-1 ;
    elseif  i == size(CN,1)
        % Last segment. We have to contemplate the following scenarios
        
        if CN(i-1,2) ~= CN(i,1) &&  CN(i,2) == 1
            % The first node of the segment is different from the last one
            % of the previous segment, and the last node of the segment is
            % equal to the first node of the list
            COORquad = [COORquad; edgePoints(1:end-1,:)];
            LINES{i} = [nodeCount + (1:(nPoints-1)),1];
        elseif  CN(i-1,2) == CN(i,1) &&  CN(i,2) == 1
            % The first node of the segment is equal to the last one
            % of the previous segment, and the last node of the segment is
            % equal to the first node of the list
            COORquad = [COORquad; edgePoints(2:end-1,:)];
            LINES{i} = [(nodeCount-1) + (1:(nPoints-1)),1];
        elseif   CN(i-1,2) ~= CN(i,1) &&  CN(i,2) ~= 1
            
            COORquad = [COORquad; edgePoints];
            LINES{i} = nodeCount + (1:nPoints);
            nodeCount = nodeCount + nPoints;
            
        else
            error('Option not contemplated')
        end
        
    end
    
    nNODES_ALL = nNODES_ALL + length( LINES{i}) ;
    
end

%


