function [PhiALL,Phi_non,G_chol] =...
    NonlinearModesNMODESgivenM(PhiMaster_lin,SNAPdisp,DOFl,ind_nonlinear,DATAoffline,G,G_chol)
%==========================================================================
% File: NonlinearModesNMODESgivenM.m
%
% Purpose
% -------
%   Build the nonlinear displacement subspace for a 1D damage manifold,
%   enforcing G-orthogonality to the linear/master mode. Given snapshots and
%   a positive-definite inner-product matrix G (e.g., geometric mass on free
%   DOFs), the routine:
%     (1) Projects nonlinear snapshots onto the G-orthogonal complement of
%         PhiMaster_lin,
%     (2) Computes a weighted (G) low-rank basis via randomized SVD,
%     (3) Truncates according to DATAoffline.errorDISP and nmodes caps,
%     (4) Returns the total basis [PhiMaster_lin, Phi_non] and reports the
%         global G-norm reconstruction error.
%
% Function Signature
% ------------------
%   [PhiALL, Phi_non, G_chol] = ...
%     NonlinearModesNMODESgivenM(PhiMaster_lin, SNAPdisp, DOFl, ind_nonlinear, ...
%                                DATAoffline, G, G_chol)
%
% Inputs
% ------
%   PhiMaster_lin : [n_free×1] linear/master mode on free DOFs (G-normalized).
%   SNAPdisp      : [n_free×Ns] displacement snapshots on free DOFs (all steps).
%   DOFl          : indices of free/independent DOFs (used to form SNAP blocks).
%   ind_nonlinear : indices of columns in SNAPdisp considered “nonlinear”.
%   DATAoffline   : struct with optional fields:
%       • errorDISP           : SVD tolerance for displacement basis (default 0).
%       • nmodes_nonlinear    : hard cap for nonlinear modes (default = rank).
%       • nmodes_INELASTIC    : final cap for inelastic modes (default = rank).
%   G             : [n_free×n_free] SPD matrix defining the inner product
%                   (e.g., geometric mass on free DOFs). REQUIRED (non-empty).
%   G_chol        : (optional) upper-triangular Cholesky of G; if empty, you
%                   should precompute it before calling for efficiency.
%
% Outputs
% -------
%   PhiALL  : [n_free×(1+r)] concatenated basis [PhiMaster_lin, Phi_non].
%   Phi_non : [n_free×r] nonlinear modes G-orthogonal to PhiMaster_lin.
%   G_chol  : Cholesky factor of G (passed-through if provided).
%
% Method (sketch)
% ---------------
%   • G-orthogonal projection:
%       SNAP_nl = SprojDEF_operator(PhiMaster_lin, G, SNAPdisp(:, ind_nonlinear))
%   • Weighted randomized SVD:
%       [U_W, S, V] = SRSVD(G_chol * SNAP_nl, DATAoffline.errorDISP)
%       Phi_non = G_chol \ U_W(:, 1:r)
%   • Mode capping:
%       r ← min(nmodes_nonlinear, nmodes_INELASTIC, length(S))
%   • Global error (diagnostic):
%       coeff = PhiALL' * G * SNAPdisp
%       errG  = ||SNAPdisp − PhiALL*coeff||_G / ||SNAPdisp||_G
%
% Notes
% -----
%   • G must be SPD on the free DOFs; set ComputeGeometricMassMatrix = true
%     in the calling input file to ensure availability.
%   • errorDISP often yields many nonlinear modes; use nmodes_nonlinear and/or
%     nmodes_INELASTIC to impose practical caps.
%   • PhiMaster_lin is assumed G-normalized and representative of the linear
%     response; SNAP selection in ind_nonlinear should exclude purely elastic
%     steps.
%
% Dependencies
% ------------
%   SprojDEF_operator, SRSVD, DefaultField
%
% Author
% ------
%   Joaquín A. Hernández Ortega (JAHO) — jhortega@cimne.upc.edu
%   CIMNE / Universitat Politècnica de Catalunya (UPC)
%
% History (versioned)
% -------------------
%   2025-11-07 — Barcelona, Spain — Header/comments generated automatically by ChatGPT; no logic changes. (JAHO)
%   2025-10-23 — Barcelona (Balmes 185), Spain — First notes on damage variant. (JAHO)
%==========================================================================

if nargin == 0
    load('tmp1.mat')
   
end

if isempty(G)
    error('You should set .... ComputeGeometricMassMatrix = true in the input data file')
end

% nonlinear COMPONENT. IT SHOULD BE ORTHOGONAL TO  PhiMaster_lin in the norm
% induced by G. The following operator computes the orthogonal projection
[SNAPdisp_nonlinear_L,~,~] =  SprojDEF_operator(PhiMaster_lin,G,SNAPdisp(DOFl,ind_nonlinear)) ;
% Basis matrix for  SNAPdisp_nonlinear_L via weighted SVD. We use randomized
% SVD for efficiency, and the norm induced by Mgeo (geometric mass matrix)
% disp(['Computing Choleski Decomposition of Geometric Mass Matrix....'])
% G_chol = chol(G) ; 
% disp(['...done'])


DATAoffline = DefaultField(DATAoffline,'errorDISP',0) ;
[Phi_non_All_W,SSS,VVV] = SRSVD(G_chol*SNAPdisp_nonlinear_L,DATAoffline.errorDISP) ;
% TRUNCATION BASED ON DATAoffline.errorDISP tends to give rather high
% number of nonlinear modes
% This is why the user should provide an additional input
% DATAoffline.nmodes_nonlinear
DATAoffline = DefaultField(DATAoffline,'nmodes_nonlinear',length(SSS)) ;
nmodes_nonlinear = min(DATAoffline.nmodes_nonlinear,length(SSS)) ;
Phi_non_All_W = Phi_non_All_W(:,1:nmodes_nonlinear) ;
Phi_non  = G_chol\Phi_non_All_W(:,1:nmodes_nonlinear) ;


DATAoffline = DefaultField(DATAoffline,'nmodes_INELASTIC',length(SSS)) ;
nmodes_INELASTIC = min(DATAoffline.nmodes_INELASTIC,length(SSS)) ;
Phi_non = Phi_non(:,1:nmodes_INELASTIC) ;
 



% DATALOCCC.Mchol = G_chol; 
% 
% Phi_non = WSVDT(Phi_non,[],DATALOCCC);
% 
%  

% Examine total error
PhiALL = [PhiMaster_lin,Phi_non] ;
coeffs = PhiALL'*G*SNAPdisp(DOFl,:) ;

ERROR_phi_mat = SNAPdisp(DOFl,:) - PhiALL*coeffs ;
nE2 = sum(sum(ERROR_phi_mat.*(G*ERROR_phi_mat))) ;
nD2 = sum(sum(SNAPdisp(DOFl,:).*(G*SNAPdisp(DOFl,:)))) ;
totalERROR = sqrt(nE2/nD2) ;
disp(['Total ERROR svd in norm G = ',num2str(totalERROR),' for ','nmodes nonlinear =',num2str(size(Phi_non,2))])


% 
% DATALOCCC.Mchol = G_chol; 
% 
% PhiALL_M = WSVDT(PhiALL,[],DATALOCCC);
% 
% coeffs = PhiALL_M'*G*SNAPdisp(DOFl,:) ;
% 
% ERROR_phi_mat = SNAPdisp(DOFl,:) - PhiALL_M*coeffs ;
% nE2 = sum(sum(ERROR_phi_mat.*(G*ERROR_phi_mat))) ;
% nD2 = sum(sum(SNAPdisp(DOFl,:).*(G*SNAPdisp(DOFl,:)))) ;
% totalERROR2 = sqrt(nE2/nD2) ;
% disp(['Total ERROR svd in norm Mll = ',num2str(totalERROR2),' for ','nmodes nonlinear =',num2str(size(Phi_non,2))])
% 
% disp(['Ration ERROR norm G /Mll = ',num2str(totalERROR/totalERROR2)])

