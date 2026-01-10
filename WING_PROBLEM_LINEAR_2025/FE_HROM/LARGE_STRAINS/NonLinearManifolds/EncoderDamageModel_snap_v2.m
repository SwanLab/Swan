function [qL_extended,qLATENT_LOC] = EncoderDamageModel_snap_v2(BasisU_master,GmatN,SNAP_cluster,DOFl,idimLAT,DATA_evaluateTAU_and_DER)
% EncoderDamageModel_snap_v2
% -------------------------------------------------------------------------
% PURPOSE
%   Encode a batch of displacement snapshots into the *v2 damage-manifold*
%   latent space and evaluate the corresponding extended coordinates τ(q).
%   Unlike the v1 encoder, **v2 does not compute qNONrel at the encoder**;
%   it simply returns q = [qLIN; qNON], and the decoder
%   `tauFUN_1paramDAMAGE_v2` safely handles qNONrel = qNON/qLIN internally.
%
% CONTEXT
%   This routine is the encoder counterpart of `tauFUN_1paramDAMAGE_v2`,
%   whose decoder is
%       τ(q) = [ qLIN ; qNON ; qLIN * g(qNON/qLIN) ],
%   with tolerance-based handling near qLIN ≈ 0. Keeping all “qLIN-division”
%   inside the decoder avoids derivative pathologies and centralizes the
%   regularization.
%
% INPUTS
%   BasisU_master  : [nDOF×n_master] master basis used for encoding,
%                    typically [Φ_lin, Φ_non,mst], G-orthonormal.
%   GmatN          : [nDOF×nDOF] SPD inner-product matrix (e.g., K_ll).
%   SNAP_cluster   : Struct holding SVD of the displacement snapshots:
%                      .DISP.U, .DISP.S, .DISP.V  so that  DISP ≈ U*S*Vᵀ.
%   DOFl           : Indices of free (independent) DOFs.
%   idimLAT        : Index of the τ component to extract from τ(q)
%                    (e.g., 1→qLIN, 2→qNON, etc.).
%   DATA_evaluateTAU_and_DER :
%                    Struct with spline evaluators and
%                    .nameFunctionEvaluate = 'tauFUN_1paramDAMAGE_v2'.
%
% OUTPUTS
%   qL_extended    : [n_modes×n_snap] extended coordinates τ(q), i.e.
%                    [qLIN ; qNON ; qLIN*g(qNON/qLIN)] stacked in modal space.
%   qLATENT_LOC    : [1×n_snap] selected component of τ(q) indexed by idimLAT.
%
% METHOD
%   1) Reconstruct snapshot matrix on the fly from its SVD factors:
%        DISP ≈ U*S*Vᵀ.
%   2) Encode by G-weighted projection onto the master basis:
%        coeffs = Φ_masterᵀ G (U*S*Vᵀ).
%        This yields q = [qLIN; qNON; ...] for the columns corresponding to
%        the master coordinates (first two rows are qLIN and qNON).
%   3) Pass q directly to the decoder:
%        qL_extended = tauFUN_1paramDAMAGE_v2(q, DATA_evaluateTAU_and_DER).
%   4) Extract requested component:
%        qLATENT_LOC = qL_extended(idimLAT,:).
%
% NUMERICAL NOTES
%   • No division by qLIN here. All handling of qNONrel = qNON/qLIN and the
%     qLIN→0 regularization is performed in the decoder, ensuring consistent
%     derivatives (τ′, τ″) across the pipeline.
%   • BasisU_master must be G-orthonormal to interpret coefficients as
%     generalized coordinates.
%
% VERSION / AUTHOR
%   24-Oct-2025 — Encoder aligned with τ_v2 (relative-coordinate handled
%   inside decoder). Keeps the encoder side simple and robust.
% -------------------------------------------------------------------------

if nargin == 0
    load('tmp1.mat')
end
coeffsMasterModes = BasisU_master'*(GmatN*SNAP_cluster.DISP.U(DOFl,:)) ;  % SVD, U,S,V
coeffsMasterModes =bsxfun(@times,coeffsMasterModes',SNAP_cluster.DISP.S)' ;
qL  = coeffsMasterModes*SNAP_cluster.DISP.V' ;
% qLIN = coeffsMasterModes(1,:) ; 
% qNON = coeffsMasterModes(2,:) ;
% 
% max_qLIN = max(abs(qLIN)) ; 
% TOL = 1e-12 ; 
% IndicesNotCloseToZero = find(abs(qLIN) > TOL*max_qLIN ) ; 
% qNONrel = zeros(size(qNON)) ; 
% qNONrel(IndicesNotCloseToZero) = qNON(IndicesNotCloseToZero)./qLIN(IndicesNotCloseToZero) ; 

%qL = [qLIN; qNON] ;
%   qL_extended = tauNON(qL) ;
[qL_extended] = feval(DATA_evaluateTAU_and_DER.nameFunctionEvaluate,qL,DATA_evaluateTAU_and_DER) ;
qLATENT_LOC = qL_extended(idimLAT,:) ;