function SNAPforceS = BasisA_InertialForces(NstRED,BasisForces,DATA)
%--------------------------------------------------------------------------
% function SNAPforceS = BasisA_InertialForces(NstRED, BasisForces, DATA)
%
% PURPOSE:
%   Computes a set of *force snapshots* corresponding to the projection 
%   of a set of reduced displacement modes (NstRED) onto a basis of 
%   internal forces (BasisForces). This operation is useful in model 
%   reduction contexts, particularly when building reduced-order models 
%   using Galerkin projection or hyperreduction strategies.
%
% INPUTS:
%   - NstRED       : [NDOF x nDEF] matrix of reduced displacement modes,
%                    where each column is a mode vector.
%   - BasisForces  : [ngausT*nstrainF x nBasisForces] matrix where each
%                    column represents a stress basis function (in reduced
%                    form), arranged per strain component and Gauss point.
%   - DATA         : Structure containing mesh information, specifically:
%       • DATA.MESH.ndim     : Number of strain components (e.g., 3 in 2D).
%       • DATA.MESH.ngausT   : Number of Gauss points in the full mesh.
%
% OUTPUT:
%   - SNAPforceS   : [ngausT x (nDEF*nBasisForces)] matrix of force 
%                    snapshots. Each block corresponds to the projection of
%                    a displacement mode over all stress basis functions.
%
% INTERNAL LOGIC:
%   For each reduced displacement mode:
%     - Loop over strain components (e.g., ε_xx, ε_yy, ε_xy in 2D).
%     - Perform a pointwise multiplication between the mode components 
%       and corresponding stress basis components.
%     - Accumulate the contributions from all strain components to 
%       compute a local force projection (SNAPloc).
%     - Append the resulting snapshot vector as a new column block in
%       SNAPforceS.
%
% REMARK:
%   The multiplication by Gauss weights (Wdom) has been removed as of 
%   12-April-2020, assuming that the BasisForces are already weighted 
%   or the weighting is handled externally.
%
% AUTHOR:
%   J.A. Hernández-Ortega, UPC, last revision on 2020-04-12
%--------------------------------------------------------------------------

if nargin == 0
    load('tmp1.mat')
end

% 

nDEF = size(NstRED,2) ; % Number of displacement modes
nBasisForces = size(BasisForces,2) ; % Number of stress modes
SNAPforceS = [] ; % Matrix of snapshots

nstrainF = DATA.MESH.ndim; 

ngausT_RHS = DATA.MESH.ngausT/DATA.MESH.ngaus_STRESS*DATA.MESH.ngaus_RHS  ;  % 30-May-2025


for I = 1:nDEF
    SNAPloc = zeros(ngausT_RHS,nBasisForces) ;
    for istrain = 1:nstrainF
        B_stress = bsxfun(@times,BasisForces(istrain:nstrainF:end,:),NstRED(istrain:nstrainF:end,I));
        % if DATAIN.NOTMULTIPLIED_BY_WEIGHTS ==0
       % SNAPloc = SNAPloc + bsxfun(@times,B_stress,sqrt(Wdom)) ;  %
       % Removed 12-APril-2020
        % else
             SNAPloc = SNAPloc + B_stress ;
        % end
    end
    SNAPforceS = [SNAPforceS SNAPloc] ;
end