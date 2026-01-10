function [inELEM,onELEM,TriLocal] = IsInsideGeneral2023(xLOC,COOR,CN,elemLOC,IND_POLYG,POLYINFO)
%--------------------------------------------------------------------------
% function [inELEM, onELEM, TriLocal] = IsInsideGeneral2023(xLOC, COOR, CN, elemLOC, IND_POLYG, POLYINFO)
%
% PURPOSE:
%   Determines whether a given point `xLOC` lies inside (or on the boundary of) 
%   a specific finite element (2D or 3D), defined by its nodal coordinates 
%   and connectivity. For 2D, a polygonal test is performed using `inpolygon`. 
%   For 3D, the function uses Delaunay triangulation and `pointLocation`.
%
% INPUTS:
%   - xLOC      : 1 x ndim array with coordinates of the query point
%   - COOR      : Nnodes x ndim array of node coordinates
%   - CN        : Nelements x Nnodes_per_elem connectivity matrix
%   - elemLOC   : Index of the element to test
%   - IND_POLYG : Indices (within the element's local nodes) of the nodes 
%                 forming the polygon/polyhedron to check against
%   - POLYINFO  : Structure with field:
%                   * TriangulationDelaunay : cell array containing cached 
%                     `delaunayTriangulation` objects for each element (optional)
%
% OUTPUTS:
%   - inELEM    : Boolean flag, 1 if point is inside the element
%   - onELEM    : Boolean flag, 1 if point lies on the boundary (2D only)
%   - TriLocal  : If 3D, local Delaunay triangulation object used for the check
%
% NOTES:
%   - In 2D: Assumes planar polygons and uses MATLAB’s `inpolygon` function.
%   - In 3D: Uses `delaunayTriangulation` and `pointLocation`. If triangulation 
%     is cached in `POLYINFO`, it is reused for efficiency.
%   - IND_POLYG is assumed to possibly include a repeated last node in 3D 
%     (for closed surfaces); this is removed prior to triangulation.
%
% DEPENDENCIES:
%   - Requires MATLAB’s Computational Geometry toolbox (for 3D check).
%
% EXAMPLE USAGE:
%   [in, on] = IsInsideGeneral2023([0.1, 0.2], COOR, CN, 5, [1 2 3 4], POLYINFO)
%
% AUTHOR:
%   Joaquín A. Hernández, UPC-CIMNE, May 2025.
%--------------------------------------------------------------------------

TriLocal = [] ; 
ndim = size(COOR,2) ;
CNloc = CN(elemLOC,:); % Nodes forming the element
if ndim == 2
    INDnodes =CNloc(IND_POLYG) ;  % Corner nodes (to define a polygon)
    COORelem = COOR(INDnodes,:) ; % Coordinates of the elemen
    [inELEM,onELEM] = inpolygon(xLOC(1),xLOC(2),COORelem(:,1),COORelem(:,2)) ;
else
    % Use triangulation
    % ------------------
    INDnodes =CNloc(IND_POLYG) ;  % Corner nodes
    INDnodes = INDnodes(1:end-1) ;
    COORelem = COOR(INDnodes,:) ; % Coordinates of the elemen
  %  error('REvise this ! ')
    if isempty(POLYINFO.TriangulationDelaunay{elemLOC})
    TriLocal = delaunayTriangulation(COORelem);
    else
        TriLocal = POLYINFO.TriangulationDelaunay{elemLOC} ; 
    end
    CHECKINSIDE = pointLocation(TriLocal,xLOC) ;
    if isnan(CHECKINSIDE)
        inELEM = 0 ; onELEM = 0 ;
    else
        inELEM = 1; onELEM = 0 ;
    end    
end