function [elemCONTAINER,POLYINFO] = WhichElementInside2023(xLOC,INDnear,VAR_SMOOTH_FE,IND_POLYG,POLYINFO,inew) 
%--------------------------------------------------------------------------
% function [elemCONTAINER, POLYINFO] = WhichElementInside2023( ...
%                     xLOC, INDnear, VAR_SMOOTH_FE, IND_POLYG, POLYINFO, inew)
%
% PURPOSE:
%   Determines which finite element (among neighbors of a given node) 
%   contains a spatial point `xLOC`. This is typically used during pointwise 
%   evaluations or mesh interpolation (e.g., in smoothing or mapping operations).
%
% INPUTS:
% -------
%   - xLOC         : [1 × ndim] vector, coordinates of the point to be located
%   - INDnear      : integer, index of the closest node to xLOC
%   - VAR_SMOOTH_FE: structure containing:
%        · CN           → Element connectivity (elements × nodes)
%        · CONNECT_info → Element neighbor information (CONNECT_info.TableElements)
%
%   - IND_POLYG    : (optional) structure for polygonal domains
%   - POLYINFO     : structure containing auxiliary triangulation data:
%        · setElements           → previous element guesses for each point
%        · TriangulationDelaunay → cell array storing local triangulations
%
%   - inew         : index of the current point (used to retrieve stored element)
%
% OUTPUTS:
% --------
%   - elemCONTAINER: index of the element that contains `xLOC`. Empty if not found.
%   - POLYINFO     : updated structure with new triangulation data (per element)
%
% METHOD:
% -------
%   1. Identify all elements that include the node `INDnear` (direct connectivity).
%
%   2. If a prior guess for the element containing `xLOC` exists in 
%      POLYINFO.setElements(inew), test that element first, along with its
%      direct neighbors (as defined in VAR_SMOOTH_FE.CONNECT_info).
%
%   3. Iterate through the list of candidate elements:
%       - Call `IsInsideALL2023()` to check if `xLOC` lies inside or on
%         the boundary of the element.
%       - If positive, set `elemCONTAINER` and break the loop.
%
%   4. Store/update local triangulation results in POLYINFO.TriangulationDelaunay.
%
% THEORETICAL CONTEXT:
% --------------------
%   This function plays a key role in **element-wise localization** needed 
%   for interpolation, projection, or integration when working with:
%     - Arbitrary geometries
%     - Pointwise quadrature selection
%     - Mesh transfer or mesh-to-mesh projection tasks
%
%   It supports adaptive cubature and empirical interpolation by helping
%   identify the finite element context for a given spatial query.
%
% EXAMPLE:
% --------
%   [elem, POLYINFO] = WhichElementInside2023(x, nodeID, VARFE, IPOLY, POLY, i);
%
% REMARKS:
% --------
%   - Uses a prioritized search strategy based on proximity and memory
%     from previous localization calls.
%   - Local triangulations are reused or updated inside POLYINFO.
%   - Depends on external function: `IsInsideALL2023`.
%
% AUTHOR:
%   Joaquín A. Hernández, UPC-CIMNE, March 2023
%--------------------------------------------------------------------------

if nargin == 0
    load('tmp.mat')
end

ielem = 1;
elemCONTAINER = [] ;
[ELEMnear, aaaa ]= find(VAR_SMOOTH_FE.CN == INDnear) ;  % Elements sharing INDnear.

if ~isempty(POLYINFO.setElements)
    currentELEMENT = POLYINFO.setElements(inew) ;
    
    if currentELEMENT ~=0
        % This is the last element stored in memory. We shall search first in
        % this element, as well as on the neighboring elements
        nNEIGHS = VAR_SMOOTH_FE.CONNECT_info.ElemShared(currentELEMENT)  ;
        NEIGH_elemes = VAR_SMOOTH_FE.CONNECT_info.TableElements(currentELEMENT,1:nNEIGHS) ;
        ELEMnear = setdiff(ELEMnear,currentELEMENT) ;
        ELEMnear = unique([ELEMnear;NEIGH_elemes(:)],'stable') ;
        % Finally
        ELEMnear = [currentELEMENT;ELEMnear] ;
        
    end
    
end

while ielem <= length(ELEMnear)
    elemLOC = ELEMnear(ielem) ;
    
    [inELEM,onELEM,TriLocal] = IsInsideALL2023(xLOC,VAR_SMOOTH_FE,elemLOC,IND_POLYG,POLYINFO)  ;
    POLYINFO.TriangulationDelaunay{elemLOC} = TriLocal ;
    if inELEM == 1 || onELEM == 1
        elemCONTAINER  = elemLOC ;
        %             if elemCONTAINER == 99
        %                 disp('Borrar esto')
        %             end
        break
    end
    ielem = ielem + 1;
end

% if isempty(elemCONTAINER)
%     disp('')
% end