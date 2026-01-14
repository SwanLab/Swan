function [COORquad, LINES] = generateInterpolatedEdges(COOR, orders)
%GENERATEINTERPOLATEDEDGES Generates nodal coordinates along multiple edges
%   [COORquad, LINES] = generateInterpolatedEdges(COOR, orders)
%
%   This function takes a set of edges, each defined by two points (start and end),
%   and generates interpolated nodal coordinates along each edge according to the 
%   specified polynomial interpolation order.
%
%   INPUT:
%     COOR   - (2N x 2) array of coordinates, where each pair of rows defines one edge:
%              row 2i-1 is the start point, row 2i is the end point for edge i.
%     orders - (N x 1) array of polynomial orders for each edge.
%              For order p, the function will place (p+1) points along the edge.
%
%   OUTPUT:
%     COORquad - (sum(p+1) x 2) array of generated coordinates for all edges.
%     LINES    - cell array where LINES{i} contains the global indices of the
%                nodes belonging to edge i, ordered from start to end.
%
%   Example usage:
%     COOR = [x1 y1; x2 y2; x3 y3; x4 y4];  % Two edges
%     orders = [4, 2];
%     [COORquad, LINES] = generateInterpolatedEdges(COOR, orders);
% JAHO + ChatGPT, 11-May-2205, Secrets by Farga

% Initialize
COORquad = [];
LINES = cell(length(orders), 1);
nodeCount = 0;

for i = 1:length(orders)
    p = orders(i);
    nPoints = p + 1;

    idx1 = 2*i - 1;  % Start point index
    idx2 = 2*i;      % End point index

    P1 = COOR(idx1, :);
    P2 = COOR(idx2, :);

    t = linspace(0, 1, nPoints)';
    edgePoints = (1 - t) * P1 + t * P2;

    COORquad = [COORquad; edgePoints];
    LINES{i} = nodeCount + (1:nPoints);
    nodeCount = nodeCount + nPoints;
end

end
