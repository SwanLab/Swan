function [wST,A] = ConvertECM_points2elements(DATA,A,wST)
function [wST,A] = ConvertECM_points2elements(DATA,A,wST)
% =========================================================================
% CONVERTECM_POINTS2ELEMENTS — Reuse point-wise ECM for element selection
% =========================================================================
% PURPOSE
%   Aggregate Gauss-point snapshots and weights into ELEMENT-wise snapshots
%   so that a point-based ECM routine can be directly reused to select
%   ELEMENTS instead of Gauss points.
%
% WHAT THIS DOES
%   • For each element e, sum the contributions of its ngausLOC Gauss points:
%       A_elem(e,:) = Σ_{g∈e}  wST(g) * A_gp(g,:)
%     where A{icluster} stores rows by Gauss-point order (element-major
%     layout g = igaus:ngaus:end). This yields one row per element.
%   • Optionally one could normalize by element “volume”
%       VolElements(e) = Σ_{g∈e} wST(g)
%     (kept here as a commented line). By default we DO NOT normalize;
%     we pass raw element aggregates to the ECM.
%   • Replace point weights by element weights for ECM. Here we set:
%       wST = ones(nelem,1)
%     so that the ECM treats all elements with equal prior weight and the
%     weight-solving step determines optimal reduced weights itself.
%
% INPUTS
%   DATA.MESH.nelem          : number of elements.
%   DATA.MESH.ngaus_STRESS   : #Gauss points per element (layout matches A).
%   A  : cell array; A{icluster} = [#(GP) × ncols] snapshots at Gauss points
%        to be used by ECM (e.g., internal-force snapshots columns).
%   wST: [#(GP) × 1] Gauss weights (include Jacobian factors if appropriate).
%
% OUTPUTS
%   wST: [nelem × 1] element weights for ECM (set to ones by design).
%   A  : cell array; A{icluster} becomes [nelem × ncols] element-aggregated
%        snapshots (sum of GP rows per element, weighted by original wST).
%
% RATIONALE / NOTES
%   • Element-based ECM is often convenient for implementation or storage
%     (fewer candidates, simpler bookkeeping of internal variables).
%   • If you require volume-averaged quantities per element, uncomment the
%     normalization line:
%         A{icluster} = bsxfun(@times, Aloc, 1./VolElements);
%     and consider setting wST = VolElements instead of ones.
%   • Ensure DATA.MESH.ngaus_STRESS matches the GP layout used to build A.
%   • The current choice wST=ones delegates all weighting to the ECM solver.
%
% COMPLEXITY
%   O(#clusters × nelem × ngausLOC) time; O(nelem) extra memory per cluster.
%
% VERSION / AUTHORSHIP
%   • 07-NOV-2025 — Commented, element-ECM adapter explained; default to
%                    unnormalized aggregation and unit element weights. Barcelona.
%   Author: Joaquín Alberto Hernández Ortega (JAHO) — jhortega@cimne.upc.edu
%   (This header was generated automatically by ChatGPT on 07-NOV-2025.)
% =========================================================================

if nargin == 0
    load('tmp.mat')
end


% Element-based ECM (elements rather than Gauss points, elements from the mesh are selected)
nelem = DATA.MESH.nelem ;
ngausLOC  = DATA.MESH.ngaus_STRESS ;
VolElements = zeros(nelem,1) ;
for igaus = 1:ngausLOC
    VolElements = VolElements + wST(igaus:ngausLOC:end) ;
end
for icluster = 1:length(A)
    Aloc = zeros(nelem,size(A{icluster},2));
    for igaus = 1:ngausLOC
        Aloc = Aloc + bsxfun(@times,A{icluster}(igaus:ngausLOC:end,:),wST(igaus:ngausLOC:end)) ;
    end
    
    %  A{icluster} = bsxfun(@times,Aloc,1./VolElements) ;
    
    A{icluster} =Aloc; % bsxfun(@times,Aloc,1./VolElements) ;
end
%wST = VolElements ;
wST = ones(size(VolElements));