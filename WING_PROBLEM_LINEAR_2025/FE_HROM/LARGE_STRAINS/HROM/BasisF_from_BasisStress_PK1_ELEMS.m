function SNAPforceS = BasisF_from_BasisStress_PK1_ELEMS(BstRED_l,BasisPone,DATA,wST)
%--------------------------------------------------------------------------
% BasisF_from_BasisStress_PK1_ELEMS
%
% PURPOSE:
%   Build the snapshot matrix of *internal force contributions* associated
%   with a given stress basis (1st Piola–Kirchhoff stresses) and reduced B
%   operators. The output collects, per element, the projection of stress
%   modes onto the force space. This is typically used in HROM/ECM contexts
%   to link stress-mode bases with reduced equilibrium equations.
%
% INPUT:
%   BstRED_l   : [nstrainF*ngauss x nDEF] reduced B operator evaluated
%                for each reduced displacement mode (columns).
%   BasisPone  : [nstrainF*ngauss x nBasisPone] basis of 1st PK stress 
%                modes (columns).
%   DATA       : struct with FE mesh information
%                  .MESH.ndim        - spatial dimension
%                  .MESH.ngausT      - total number of Gauss points
%                  .MESH.ngaus_STRESS- number of Gauss points per element
%                  .MESH.nelem       - number of elements
%   wST        : [ngauss x 1] weights associated with the stress Gauss points.
%
% OUTPUT:
%   SNAPforceS : [nelem x (nDEF*nBasisPone)] matrix of elemental internal
%                force snapshots, each column corresponding to the
%                force pattern induced by one (stress-mode × displacement-mode) 
%                combination.
%
% METHOD:
%   1. For each displacement mode, loop over all stress components.
%   2. Multiply stress basis functions by the reduced B-operator entries
%      (BstRED_l) and Gauss weights (wST), accumulating Gauss-point
%      contributions.
%   3. Assemble snapshots per element by summing Gauss-point contributions
%      within each element.
%   4. Concatenate results across displacement modes.
%
% NOTES:
%   - nstrainF = ndim^2 (Voigt-like indexing for 1st PK stress tensor).
%   - Output is organized per element, not per Gauss point.
%   - Provides the stress-to-force mapping needed for reduced internal force
%     evaluation in hyperreduction.
%--------------------------------------------------------------------------

if nargin == 0
    load('tmp1.mat')
end

% 

nDEF = size(BstRED_l,2) ; % Number of displacement modes
nBasisPone = size(BasisPone,2) ; % Number of stress modes
SNAPforceS = [] ; % Matrix of snapshots
nstrainF = DATA.MESH.ndim^2; 

for I = 1:nDEF
    SNAPloc = zeros(DATA.MESH.ngausT,nBasisPone) ;
    for istrain = 1:nstrainF
        B_stress = bsxfun(@times,BasisPone(istrain:nstrainF:end,:),BstRED_l(istrain:nstrainF:end,I));
        
             SNAPloc = SNAPloc + bsxfun(@times,B_stress,wST) ; 
        
    end
    SNAPforceS = [SNAPforceS SNAPloc] ;
end


% Summing up the contribution of each element 
nelem = DATA.MESH.nelem ; 

SNAPforceS_elem = zeros(nelem,size(SNAPforceS,2));
ngausLOC  = DATA.MESH.ngaus_STRESS ; 

for igaus = 1:ngausLOC
    SNAPforceS_elem = SNAPforceS_elem + SNAPforceS(igaus:ngausLOC:end,:) ; 
end


SNAPforceS = SNAPforceS_elem  ; 
