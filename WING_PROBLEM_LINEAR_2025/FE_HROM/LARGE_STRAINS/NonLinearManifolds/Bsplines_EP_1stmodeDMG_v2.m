function [PhiNONmst,PhiNONslv,qNON,...
    gNONslv,...
    DATA_evaluateTAU_and_DER] = Bsplines_EP_1stmodeDMG_v2(PhiLIN,Phi_non,GmatN,DATA_interp,SNAPdisp_L,GmatN_chol)
% Bsplines_EP_1stmodeDMG_v2
% -------------------------------------------------------------------------
% PURPOSE
%   Damage-aware variant of the “first-mode master” manifold builder that
%   prepares the *encoder/decoder* consistent with the v2 τ-map:
%
%       τ(q) = [ qLIN ;
%                 qNON ;
%                 qLIN * g(qNONrel) ],    with  qNONrel = qNON / qLIN  (set qNONrel=0 if qLIN≈0).
%
%   From linear/nonlinear modal snapshots it returns:
%     (i)   the nonlinear *master* mode Φ_non,mst and its *slave* complement Φ_non,slv,
%     (ii)  the master coordinate qNON and the normalized slave amplitudes gNONslv,
%     (iii) spline evaluators for g(·), g′(·), g″(·) as functions of qNONrel.
%
% WHAT CHANGED VS. v1
%   • The regression abscissa is the *relative* coordinate qNONrel = qNON/qLIN
%     (instead of qNON). This makes τ and its derivatives well-posed near qLIN=0
%     and matches the v2 decoder used in tauFUN_1paramDAMAGE_v2.
%
% WORKFLOW
%   1) Linear coordinate:
%        qLIN = Φ_linᵀ G SNAP_L.  (G := GmatN; Φ_lin is G-orthonormal)
%
%   2) Nonlinear raw coordinates and master selection:
%        Q_non_all = Φ_nonᵀ G SNAP_L.
%        Normalize by qLIN:  Q_non_all./qLIN  (columnwise) to remove scaling.
%        Choose a *monotone* linear combination of the normalized rows
%        (ChooseMonotonicLinearCombination) → define ΦNONmst; G-normalize it.
%
%   3) Master coordinate and relative abscissa:
%        qNON    = ΦNONmstᵀ G SNAP_L.
%        qNONrel = qNON ./ qLIN  (elementwise); if |qLIN|<TOL, set qNONrel=0.
%
%   4) Slave basis and amplitudes:
%        Build G-orthogonal complement of ΦNONmst in Φ_non → ΦNONslv (re-orthonormalize with WSVDT).
%        qNONslv = ΦNONslvᵀ G SNAP_L.
%        gNONslv = qNONslv ./ qLIN.        % matches τ’s factorization qLIN*g(·)
%
%   5) B-spline learning in compact SVD space:
%        Perform SVDT on gNONslv to obtain U,S,V (compact slave coordinates).
%        Fit each row of Vᵀ as a function of qNONrel with cubic B-splines,
%        optionally after uniform subsampling in qNONrel for robust knot placement.
%        Package spline handles/metadata in DATA_evaluateTAU_and_DER
%        so tauFUN_1paramDAMAGE_v2 can evaluate g, g′, g″ consistently.
%
% INPUTS
%   PhiLIN      : [nDOF×1] linear/elastic mode (G-orthonormal).
%   Phi_non     : [nDOF×r_nl] nonlinear basis (G-orthonormal columns).
%   GmatN       : [nDOF×nDOF] SPD inner-product matrix (e.g., K_ll).
%   DATA_interp : Options for spline fitting and subsampling:
%                 .NSAMPLES, .order_Bslines, .ratio_NSAMPLES_knots,
%                 .SubSampling_qINELASTIC_master (logical), etc.
%   SNAPdisp_L  : [nDOF×n_snap] snapshot matrix on free DOFs.
%   GmatN_chol  : Cholesky-like factor of GmatN for weighted SVD routines.
%
% OUTPUTS
%   PhiNONmst   : [nDOF×1] nonlinear *master* mode.
%   PhiNONslv   : [nDOF×(r_nl−1)] G-orthonormal *slave* modes.
%   qNON        : [1×n_snap] master coordinate.
%   gNONslv     : [(r_nl−1)×n_snap] slave amplitudes normalized by qLIN.
%   DATA_evaluateTAU_and_DER :
%                 Struct with spline evaluators (g, g′, g″) vs qNONrel and
%                 metadata (bounds, extrapolation) for τ, τ′, τ″ assembly.
%
% NUMERICAL NOTES
%   • Division by qLIN is safeguarded: when |qLIN|<TOL, use qNONrel=0. This
%     matches the true behavior τ_i(qLIN=0,·)=0 and avoids spurious blow-ups
%     in derivatives.
%   • Using gNONslv = qNONslv./qLIN aligns training targets with the v2 decoder.
%   • Working in the compact SVD of gNONslv minimizes the number of splines
%     to fit while preserving accuracy.
%
% DELIVERABLE
%   Manifold ingredients (Φ_non,mst/slv, qNON, gNONslv, spline handles) ready
%   for tauFUN_1paramDAMAGE_v2 and for consistent-tangent assembly in ROM/HROM.
%
% VERSION / AUTHOR
%   24-Oct-2025 — J.A. Hernández. v2 encoder/decoder alignment using qNONrel.
% Comments by ChatGPT-5
% -------------------------------------------------------------------------


if nargin == 0
    load('tmp.mat')
    %       DATA_interp.NSAMPLES = 300 ;
    %       DATA_interp.ratio_NSAMPLES_knots = 0.8;
    %       DATA_interp.SubSampling_qnonlinear_master = true;
    close all
end

%
% % First latent coordinate
qLIN = PhiLIN'*GmatN*SNAPdisp_L ;
figure(190)
hold on
xlabel('Step in inelastic loading')
ylabel('Linear latent variable')
plot(qLIN,'DisplayName','Original')

% Remaining nonlinear modes, divided by qLIN
% -----------------------------------------
qNONall =  Phi_non'*GmatN*SNAPdisp_L ;
qNONall_q1 =  bsxfun(@times,qNONall',1./qLIN')' ;

%  Linear combination of qNONall_q1 which is more "monotonous"
[~,coeffsLINEAR_COMBINATION] = ChooseMonotonicLinearCombination(qNONall_q1') ;


% Accordingly, the nonlinear master mode is:
PhiNONmst = Phi_non*coeffsLINEAR_COMBINATION ;
PhiNONmst_norm = sqrt(PhiNONmst'*GmatN*PhiNONmst) ;
PhiNONmst = PhiNONmst/PhiNONmst_norm ;

% Accordintly, the nonlinear latent variable is:
qNON = PhiNONmst'*GmatN*SNAPdisp_L ;
%[qNON,IndicesNonRepeated,bbb ]= unique(qNON) ; % REmove repeated indices
%qLIN = qLIN(IndicesNonRepeated) ;
% On the other hand,
% \qNONrel = \dfrac{\qNON}{\qLIN}
qNONrel  =  bsxfun(@times,qNON',1./qLIN')' ;

signqNONrel_training = unique(sign(qNONrel)) ; 

if length(signqNONrel_training) > 1 
    error('The same should be the same for all samples ')
end

%
figure(191)
title('Nonlinear latent variable (qNON), and qNONrel as a function of the loading parameter')
hold on
xlabel('Step in inelastic loading')
ylabel('Modal amplitude   nonlinear latent variable')

plot(qNON,'DisplayName','qNON' )
plot(qNONrel,'DisplayName','qNONrel = qNON/qLIN' )

legend show

% Next step: orthogonal complement of  PhiNONmst
[PhiNONslv,~,~] =  SprojDEF_operator(PhiNONmst,GmatN,Phi_non) ;
DATAlocc.Mchol = GmatN_chol ;
DATAlocc.RELATIVE_SVD = 1;
DATAlocc.TOL =1e-10 ;
[PhiNONslv,SSs,VVv] = WSVDT(PhiNONslv,[],DATAlocc) ;

% Thus
qNONslv = PhiNONslv'*GmatN*SNAPdisp_L  ;
% Next we divided qNONslv by qLIN
gNONslv =  bsxfun(@times,qNONslv',1./qLIN')' ;

SHOW_before_SVD = 1;

if SHOW_before_SVD == 1
    figure(456)
    hold on
    title(['Slave modes amplituted divided by qLIN (gNONslv) as a function of qNONrel = qNON/qLIN '])
    for iii  =1:size(gNONslv,1)
        plot(qNONrel,gNONslv(iii,:),'DisplayName',['gNONslv_{',num2str(iii),'}'])
    end
    legend show
end



DATA_interp = DefaultField(DATA_interp,'SubSampling_qINELASTIC_master',true) ;

if ~DATA_interp.SubSampling_qINELASTIC_master
    
    [UU,SS,VV] = SVDT(gNONslv) ;
    [DATA_evaluateTAU_and_DER] = BsplinesLeastSquares_DAMAGE(DATA_interp, qNONrel, VV', UU, SS) ;
    
else
    [DATA_evaluateTAU_and_DER] =  BsplinesSUBSAMPL_damage(qNONrel,DATA_interp,gNONslv) ;
    
    
    
end
DATA_evaluateTAU_and_DER.signqNONrel_training = signqNONrel_training ; 

