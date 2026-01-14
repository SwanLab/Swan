function [inELEM,onELEM,TriLocal] = IsInsideALL2023(xLOC,VAR_SMOOTH_FE,elemLOC,IND_POLYG,POLYINFO)
%--------------------------------------------------------------------------
% function [inELEM, onELEM, TriLocal] = IsInsideALL2023(xLOC, VAR_SMOOTH_FE, elemLOC, IND_POLYG, POLYINFO)
%
% PURPOSE:
%   Unified function to determine whether a point `xLOC` lies inside a given 
%   finite element `elemLOC`, accounting for 1D, 2D, or 3D geometries.
%   In 1D, it performs a simple interval check.
%   In 2D/3D, it delegates the check to `IsInsideGeneral2023`.
%
% INPUTS:
%   - xLOC       : Coordinates of the query point (1x1 for 1D, 1x2 for 2D, 1x3 for 3D)
%   - VAR_SMOOTH_FE : Structure with fields:
%                     * COOR: coordinates of all nodes
%                     * CN  : connectivity matrix (Nelements x Nnodes_per_elem)
%   - elemLOC    : Index of the element to test
%   - IND_POLYG  : Indices of corner nodes (within local connectivity) to define polygon/polyhedron
%   - POLYINFO   : Structure containing optional cached triangulations for 3D case
%
% OUTPUTS:
%   - inELEM     : Boolean flag, 1 if the point is inside the element
%   - onELEM     : Boolean flag, 1 if the point lies on the boundary (2D only)
%   - TriLocal   : Local triangulation object used (non-empty only in 3D)
%
% NOTES:
%   - For 1D: assumes linear elements with 2 nodes per element.
%   - For ndim > 1: uses helper function `IsInsideGeneral2023`.
%
% EXAMPLE USAGE:
%   [in, on] = IsInsideALL2023([0.3 0.2], VAR_SMOOTH_FE, 7, [1 2 3], POLYINFO);
%
% AUTHOR:
%   Joaquín A. Hernández, UPC-CIMNE, May 2025.
%--------------------------------------------------------------------------

TriLocal = [] ; 
if length(xLOC) ==1
    % 1D problem
    elemCN = VAR_SMOOTH_FE.CN(elemLOC,:) ;
    COORelemLOC = VAR_SMOOTH_FE.COOR(elemCN,:) ;
    
    inELEM = 0 ; onELEM = 0 ;
    if  xLOC >= COORelemLOC(1) && xLOC <= COORelemLOC(2)
        inELEM = 1;
    end
else
    [inELEM,onELEM,TriLocal] = IsInsideGeneral2023(xLOC,VAR_SMOOTH_FE.COOR,VAR_SMOOTH_FE.CN,elemLOC,IND_POLYG,POLYINFO) ;
end