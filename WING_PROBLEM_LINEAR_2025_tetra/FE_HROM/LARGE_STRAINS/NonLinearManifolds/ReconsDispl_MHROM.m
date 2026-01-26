function [DISP_LOC_HROM,qL] = ReconsDispl_MHROM(DISP_LOC_HROM,DISP_CONDITIONS_HROM,DATA_evaluateTAU_and_DER,BasisU,DISP_CONDITIONS,...
    BasisU_r)
% ReconsDispl_MHROM
%
% PURPOSE
%   Reconstruct **full-order nodal displacements** from the manifold HROM
%   state by:
%     1) splitting the ROM coordinates into free (qL) and boundary blocks (qR),
%     2) **extending** the free block through the encoder τ(q):  qL → τ(qL),
%     3) assembling the physical displacement with the **affine BC** basis
%        shown in the slides (screenshot / Eq. (257)):
%
%          Φ̄ = [AΦ  0] + [0 0; 0 D^U]  ⇒  d = (AΦ)·τ(qL) + D^U·qR .
%
%   This enforces essential (Dirichlet) BCs in the trial space, consistent with
%   the weak-form derivation and free/prescribed block split (Kll/Klr logic). 
%
% INPUTS
%   DISP_LOC_HROM          : Current ROM state [ qL ; qR ] (before extension).
%   DISP_CONDITIONS_HROM   : Struct with fields .DOFl, .DOFr (ROM partition).
%   DATA_evaluateTAU_and_DER: Precompiled τ evaluator (from B-spline LS), used as
%                             feval(NAME, qL, DATA) → τ(qL). 
%   BasisU                 : Φ (free-DOF basis).
%   DISP_CONDITIONS        : Must contain .A (affine inserter) and .DOFr (FE ids).
%   BasisU_r               : D^U (Dirichlet/boundary columns on DOFr).
%
% OUTPUTS
%   DISP_LOC_HROM          : Full-order displacement(s), assembled as:
%                              d = A*Φ*τ(qL) + D^U*qR  (affine case), or
%                              d(DOFl)=Φ*τ(qL), d(DOFr)=D^U*qR (standard case).
%   qL                     : Returned for convenience (free, pre-extended coords).
%
% THEORY TIE-INS (slides)
%   • Essential vs natural BCs: trial functions satisfy u=g at Dirichlet nodes; the
%     unknowns are the free block (l). :contentReference[oaicite:1]{index=1}
%   • Block elimination form Kll dl = Fl − Klr dr motivates the [qL;qR] split. :contentReference[oaicite:2]{index=2}
%   • Vectorized operators (B, L, N, ω) consume the **full** d reconstructed here. :contentReference[oaicite:3]{index=3}
%
% NOTES
%   • If DISP_CONDITIONS.A exists, we assemble with the **affine** map A*Φ and
%     add the boundary component D^U*qR on DOFr; otherwise we place Φ and D^U
%     directly on DOFl/DOFr.
%   • τ(q) may be nonlinear; this routine only needs τ (not τ′) for reconstruction.
% JAHO, 8-Oct-2025, Wednesday, 20:21, Starbucks Paseo de Gracia, BArcelona
% Comments by ChatGPT5
%==========================================================================

if nargin ==0
    load('tmp1.mat')
end

% Reconstruction displacements
% Extended coordinates
qL = DISP_LOC_HROM(DISP_CONDITIONS_HROM.DOFl,:) ;
qR = DISP_LOC_HROM(DISP_CONDITIONS_HROM.DOFr,:) ;
%   qL_extend = tauNON(qL) ;
qL_extend = feval(DATA_evaluateTAU_and_DER.nameFunctionEvaluate,qL,DATA_evaluateTAU_and_DER) ;
DISP_LOC_HROM_proj= [qL_extend;qR] ;

DOFlHROM = 1:size(BasisU,2) ; DOFrHROM = (size(BasisU,2)+1):size(DISP_LOC_HROM_proj,1);
ndof = length(DISP_CONDITIONS.DOFl)+ length(DISP_CONDITIONS.DOFr);
DISP_LOC_HROM = zeros(ndof,size(DISP_LOC_HROM_proj,2));
if  isfield(DISP_CONDITIONS,'A')
    % Constructing BasisUall
    
    
    DISP_LOC_HROM  = DISP_CONDITIONS.A*BasisU*DISP_LOC_HROM_proj(DOFlHROM,:) ;
    DISP_LOC_HROM(DISP_CONDITIONS.DOFr,:) =DISP_LOC_HROM(DISP_CONDITIONS.DOFr,:) + BasisU_r*DISP_LOC_HROM_proj(DOFrHROM,:) ;
    
    
else
    
    
    % Constructing BasisUall
    % DOFlHROM = 1:size(BasisU,2) ; DOFrHROM = (size(BasisU,2)+1):size(DISP_LOC_HROM_proj,1);
    % ndof = length(DISP_CONDITIONS.DOFl)+ length(DISP_CONDITIONS.DOFr);
    % DISP_LOC_HROM = zeros(ndof,size(DISP_LOC_HROM_proj,2));
    
    DISP_LOC_HROM(DISP_CONDITIONS.DOFl,:) = BasisU*DISP_LOC_HROM_proj(DOFlHROM,:) ;
    DISP_LOC_HROM(DISP_CONDITIONS.DOFr,:) = BasisU_r*DISP_LOC_HROM_proj(DOFrHROM,:) ;
    
end