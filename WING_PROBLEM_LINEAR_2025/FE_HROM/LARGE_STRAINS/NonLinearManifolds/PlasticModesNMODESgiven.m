function [PhiALL,Phi_non,SNAPdisp_plast_L] =...
    PlasticModesNMODESgiven(PhiMaster_lin,SNAPdisp,DOFl,ind_plastic,DATAoffline,Kll_chol,Kll) 
% =========================================================================
% PLASTICMODESNMODESGIVEN — Kll-weighted plastic basis via randomized SVD
% =========================================================================
% PURPOSE
%   From displacement snapshots restricted to PLASTIC states, build a plastic
%   basis that is Kll-orthogonal to the elastic mode(s). Steps:
%     1) Kll-orthogonal projection of plastic snapshots onto the complement of
%        PhiMaster_lin using SprojDEF_operator.
%     2) Weighted randomized SVD (SRSVD) of Kll_chol * SNAPdisp_plast_L for
%        efficiency and numerical robustness.
%     3) Truncate by DATAoffline.nmodes_PLASTIC (cap), in addition to the
%        tolerance DATAoffline.errorDISP used inside SRSVD.
%     4) Map weighted left singular vectors back to displacement space via
%        Kll_chol\· to obtain the plastic basis Φ_non.
%     5) Report global Kll-norm reconstruction error with Φ_ALL = [Φ_el, Φ_non].
%
% INPUTS
%   PhiMaster_lin : Elastic basis on DOFl (columns Kll-orthonormal).
%   SNAPdisp      : [nDOF × nsnap] displacement snapshots (all states).
%   DOFl          : Indices of free/unconstrained DOFs.
%   ind_plastic   : Column indices selecting PLASTIC snapshots.
%   DATAoffline   : Struct with fields:
%                     • errorDISP        – truncation tol for SRSVD (default 0).
%                     • nmodes_PLASTIC   – cap on # plastic modes (default = rank).
%   Kll_chol      : Cholesky-like factor s.t. Kll ≈ Kll_chol' * Kll_chol.
%   Kll           : Constrained stiffness on DOFl (inner-product metric).
%
% OUTPUTS
%   PhiALL           : [Φ_elastic, Φ_non] concatenated basis.
%   Phi_non          : Plastic basis (columns Kll-orthonormal to Φ_elastic).
%   SNAPdisp_plast_L : Kll-projected plastic snapshots used in SVD.
%
% METHOD
%   • Projection to plastic subspace:
%       [SNAPdisp_plast_L,~,~] = SprojDEF_operator(PhiMaster_lin, Kll, SNAPdisp(DOFl, ind_plastic));
%     which removes the Kll-component along Φ_elastic.
%   • Weighted randomized SVD:
%       [U_W, S, V] = SRSVD(Kll_chol * SNAPdisp_plast_L, DATAoffline.errorDISP);
%     Truncate to nmodes_PLASTIC, then unweight: Φ_non = Kll_chol \ U_W(:,1:m).
%   • Quality check (global):
%       coeff = Φ_ALL' * Kll * SNAPdisp(DOFl,:);
%       err   = || SNAPdisp(DOFl,:) − Φ_ALL*coeff ||_Kll / || SNAPdisp(DOFl,:) ||_Kll
%       where ||X||_Kll^2 = trace(X' * Kll * X).
%
% DIAGNOSTICS / NOTES
%   • DATAoffline.errorDISP alone can retain many modes; use nmodes_PLASTIC to
%     cap dimensionality explicitly.
%   • The Kll weighting (via Kll_chol) is essential: performing SVD in the
%     energy norm aligns the basis with stress/strain accuracy criteria.
%   • Outputs SingValuesPlast (local variable) are normalized by S(1) for
%     quick inspection (not returned).
%
% DEPENDENCIES
%   SprojDEF_operator, SRSVD, DefaultField.
%
% VERSION / AUTHORSHIP
%   17-AUG-2025 — Initial structured version aligned with Kll-weighted SVD,
%                 Molinos Marfagones (Cartagena).
%   07-NOV-2025 — Comments refreshed/clarified; projection/SVD pipeline and
%                 error metric spelled out. Barcelona.
%   Author: Joaquín Alberto Hernández Ortega (JAHO) — jhortega@cimne.upc.edu
%   (This header was generated automatically by ChatGPT on 07-NOV-2025.)
% =========================================================================

if nargin == 0
    load('tmp1.mat')
    DATAoffline.nmodes_PLASTIC = 50 ;
    
end

% PLASTIC COMPONENT. IT SHOULD BE ORTHOGONAL TO  PhiMaster_lin in the norm
% induced by Kll. The following operator computes the orthogonal projection
[SNAPdisp_plast_L,~,~] =  SprojDEF_operator(PhiMaster_lin,Kll,SNAPdisp(DOFl,ind_plastic)) ;
% Basis matrix for  SNAPdisp_plast_L via weighted SVD. We use randomized
% SVD for efficiency
DATAoffline = DefaultField(DATAoffline,'errorDISP',0) ;
[Phi_non_All_W,SSS,VVV] = SRSVD(Kll_chol*SNAPdisp_plast_L,DATAoffline.errorDISP) ;
% TRUNCATION BASED ON DATAoffline.errorDISP tends to give rather high
% number of plastic modes
% This is why the user should provide an additional input
% DATAoffline.nmodes_PLASTIC
DATAoffline = DefaultField(DATAoffline,'nmodes_PLASTIC',length(SSS)) ;
nmodes_PLASTIC = min(DATAoffline.nmodes_PLASTIC,length(SSS)) ;
Phi_non_All_W = Phi_non_All_W(:,1:nmodes_PLASTIC) ;
Phi_non  = Kll_chol\Phi_non_All_W(:,1:nmodes_PLASTIC) ;
SingValuesPlast = SSS(1:nmodes_PLASTIC)/SSS(1) ;

% Examine total error
PhiALL = [PhiMaster_lin,Phi_non] ;
coeffs = PhiALL'*Kll*SNAPdisp(DOFl,:) ;

ERROR_phi_mat = SNAPdisp(DOFl,:) - PhiALL*coeffs ;
nE2 = sum(sum(ERROR_phi_mat.*(Kll*ERROR_phi_mat))) ;
nD2 = sum(sum(SNAPdisp(DOFl,:).*(Kll*SNAPdisp(DOFl,:)))) ;
totalERROR = sqrt(nE2/nD2) ;
disp(['Total ERROR svd in norm Kll = ',num2str(totalERROR),' for ','nmodes PLASTIC =',num2str(size(Phi_non,2))])


