function [PhiALL,Phi_non,SNAPdisp_plast_L] =...
    PlasticModesNMODESgivenM(PhiMaster_lin,SNAPdisp,DOFl,ind_plastic,DATAoffline,Kll_chol,Kll,MgeoLL)
% =========================================================================
% PLASTICMODESNMODESGIVENM — Plastic basis with geometric-mass (Mgeo) metric
% =========================================================================
% PURPOSE
%   Build a PLASTIC reduced basis from plastic snapshots using a geometric
%   mass–weighted SVD, while enforcing Kll-orthogonality to the elastic mode.
%   Steps:
%     1) Kll-orthogonal projection of plastic snapshots off Φ_elastic.
%     2) Randomized SVD on Mgeo-weighted snapshots (chol(MgeoLL)*·) for
%        efficient, metric-aware compression.
%     3) Map weighted vectors back to displacement space and re-orthonormalize
%        in the Kll metric (WSVDT) to ensure energy-orthogonality.
%     4) Report global reconstruction errors in both Kll and Mgeo norms.
%
% CONTEXT / RATIONALE
%   • Some applications favour the geometric mass metric Mgeo for plastic
%     content extraction (e.g., force-driven plastic response), but online
%     residuals and stresses are evaluated in the Kll (stiffness) metric.
%     Hence the two-stage process: Mgeo-SVD for extraction + Kll re-orthonormalization.
%
% INPUTS
%   PhiMaster_lin : Elastic basis on DOFl (Kll-orthonormal columns).
%   SNAPdisp      : [nDOF × nsnap] displacement snapshots (all states).
%   DOFl          : Indices of free/unconstrained DOFs.
%   ind_plastic   : Column indices selecting PLASTIC snapshots.
%   DATAoffline   : Struct with:
%                     • errorDISP        — SRSVD truncation tolerance (default 0)
%                     • nmodes_PLASTIC   — hard cap on # plastic modes
%   Kll_chol      : Cholesky-like factor with Kll ≈ Kll_chol' * Kll_chol.
%   Kll           : Constrained stiffness (energy metric) on DOFl.
%   MgeoLL        : Geometric mass matrix on DOFl (SPD). REQUIRED.
%
% OUTPUTS
%   PhiALL           : [Φ_elastic, Φ_plastic] concatenated basis (Kll-orthonormal).
%   Phi_non          : Plastic basis (columns re-orthonormalized in Kll).
%   SNAPdisp_plast_L : Kll-projected plastic snapshots used for SVD.
%
% METHOD (pipeline)
%   1) Projection in Kll:
%        [SNAPdisp_plast_L,~,~] = SprojDEF_operator(PhiMaster_lin, Kll, SNAPdisp(DOFl, ind_plastic));
%   2) Mgeo-weighted randomized SVD:
%        MgeoLL_chol = chol(MgeoLL);
%        [U_W, S, V] = SRSVD(MgeoLL_chol * SNAPdisp_plast_L, DATAoffline.errorDISP);
%        Truncate to m = min(nmodes_PLASTIC, rank) and unweight:
%        Φ_non_raw = MgeoLL_chol \ U_W(:,1:m);
%   3) Re-orthonormalize in Kll (energy metric) for online consistency:
%        Φ_non = WSVDT(Φ_non_raw, [], struct('Mchol', Kll_chol));
%   4) Diagnostics:
%        • Compute global reconstruction error in Kll with Φ_ALL = [Φ_el, Φ_non].
%        • Repeat in Mgeo (re-orthonormalize Φ_ALL with MgeoLL_chol via WSVDT)
%          and print both errors and their ratio.
%
% PRACTICAL NOTES / GOTCHAS
%   • MgeoLL must be SPD and consistent with DOFl; set
%       ComputeGeometricMassMatrix = true
%     in the input pipeline producing OPERFE.Mgeo.
%   • DATAoffline.errorDISP alone may retain many modes; use nmodes_PLASTIC
%     to cap the plastic subspace dimension explicitly.
%   • Final basis returned is Kll-orthonormal (suitable for energy-norm
%     projections and online residual assembly).
%
% DEPENDENCIES
%   SprojDEF_operator, SRSVD, WSVDT, DefaultField, chol.
%
% VERSION HISTORY / AUTHORSHIP
%   • 18-OCT-2025 — Initial Mgeo-weighted variant. Saturday, 1 Balmes 185, Barcelona.
%   • 07-NOV-2025 — Comments refreshed; clarified dual-metric workflow and
%                    error reporting in Kll vs Mgeo. Barcelona.
%   Author: Joaquín Alberto Hernández Ortega (JAHO) — jhortega@cimne.upc.edu
%   (This header was generated automatically by ChatGPT on 07-NOV-2025.)
% =========================================================================

% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/112_NonLIN_ROM_RBF/15_FORCE_plast.mlx
% JAHO, 18-Oct-2025, Saturday, 1 Balmes 185, Barcelona
if nargin == 0
    load('tmp1.mat')
    DATAoffline.nmodes_PLASTIC = 1e20 ; 
    DATAoffline.errorDISP = 0 ; 
end

if isempty(MgeoLL)
    error('You should set .... ComputeGeometricMassMatrix = true in the input data file')
end

% PLASTIC COMPONENT. IT SHOULD BE ORTHOGONAL TO  PhiMaster_lin in the norm
% induced by Kll. The following operator computes the orthogonal projection
[SNAPdisp_plast_L,~,~] =  SprojDEF_operator(PhiMaster_lin,Kll,SNAPdisp(DOFl,ind_plastic)) ;
% Basis matrix for  SNAPdisp_plast_L via weighted SVD. We use randomized
% SVD for efficiency, and the norm induced by Mgeo (geometric mass matrix)
disp(['Computing Choleski Decomposition of Geometric Mass Matrix....'])
MgeoLL_chol = chol(MgeoLL) ; 
disp(['...done'])


DATAoffline = DefaultField(DATAoffline,'errorDISP',0) ;
[Phi_non_All_W,SSS,VVV] = SRSVD(MgeoLL_chol*SNAPdisp_plast_L,DATAoffline.errorDISP) ;
% TRUNCATION BASED ON DATAoffline.errorDISP tends to give rather high
% number of plastic modes
% This is why the user should provide an additional input
% DATAoffline.nmodes_PLASTIC
DATAoffline = DefaultField(DATAoffline,'nmodes_PLASTIC',length(SSS)) ;
nmodes_PLASTIC = min(DATAoffline.nmodes_PLASTIC,length(SSS)) ;
Phi_non_All_W = Phi_non_All_W(:,1:nmodes_PLASTIC) ;
Phi_non  = MgeoLL_chol\Phi_non_All_W(:,1:nmodes_PLASTIC) ;

DATALOCCC.Mchol = Kll_chol; 

Phi_non = WSVDT(Phi_non,[],DATALOCCC);

 

% Examine total error
PhiALL = [PhiMaster_lin,Phi_non] ;
coeffs = PhiALL'*Kll*SNAPdisp(DOFl,:) ;

ERROR_phi_mat = SNAPdisp(DOFl,:) - PhiALL*coeffs ;
nE2 = sum(sum(ERROR_phi_mat.*(Kll*ERROR_phi_mat))) ;
nD2 = sum(sum(SNAPdisp(DOFl,:).*(Kll*SNAPdisp(DOFl,:)))) ;
totalERROR = sqrt(nE2/nD2) ;
disp(['Total ERROR svd in norm Kll = ',num2str(totalERROR),' for ','nmodes PLASTIC =',num2str(size(Phi_non,2))])



DATALOCCC.Mchol = MgeoLL_chol; 

PhiALL_M = WSVDT(PhiALL,[],DATALOCCC);

coeffs = PhiALL_M'*MgeoLL*SNAPdisp(DOFl,:) ;

ERROR_phi_mat = SNAPdisp(DOFl,:) - PhiALL_M*coeffs ;
nE2 = sum(sum(ERROR_phi_mat.*(MgeoLL*ERROR_phi_mat))) ;
nD2 = sum(sum(SNAPdisp(DOFl,:).*(MgeoLL*SNAPdisp(DOFl,:)))) ;
totalERROR2 = sqrt(nE2/nD2) ;
disp(['Total ERROR svd in norm Mll = ',num2str(totalERROR2),' for ','nmodes PLASTIC =',num2str(size(Phi_non,2))])

disp(['Ration ERROR norm Kll /Mll = ',num2str(totalERROR/totalERROR2)])

