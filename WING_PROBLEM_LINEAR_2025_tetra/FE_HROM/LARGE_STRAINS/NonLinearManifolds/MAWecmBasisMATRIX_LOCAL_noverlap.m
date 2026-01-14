function [U,b,wADAPT,DATAoffline,numberMODES] = MAWecmBasisMATRIX_LOCAL_noverlap(Ufixed,zINI,qLATENT,A,DATAoffline,wINI,Wfe,Q)
%==========================================================================
% MAWecmBasisMATRIX_LOCAL_noverlap
%
% PURPOSE
% -------
% Build, for each latent cluster, the local exactness basis U{c} (restricted
% to candidate ECM points zINI), the target integrals b{c}, and initialize
% the per-cluster adaptive weights wADAPT(:,c)=wINI. This is the **no-overlap**
% variant: each cluster is treated independently (no shared samples).
%
% CONTEXT
% -------
% Used by MAW-ECM / CECM “PLaCE” pipelines to (i) inject a fixed elastic
% scaffold Ufixed (elastic modes + constant function), (ii) optionally
% denoise plastic/damage snapshots through a Wfe-orthonormal Q, and (iii)
% produce restricted bases and right-hand sides consistent with the current
% ECM candidate set zINI.
%
% INPUTS
% ------
% Ufixed     : Fixed scaffold common to all clusters (e.g., SVDT([Uel, ones])).
% zINI       : Vector of candidate quadrature-point indices (global → local).
% qLATENT    : Latent coordinates (only its length is used here for sizing).
% A          : 1×nC cell; A{c} ∈ R^(nIP×m_c), snapshots for cluster c.
% DATAoffline: Struct with options. Recognized field:
%              • Project_SNAP_withQ__MAW_ECM (default 1):
%                   1 → project snapshots with Q(QᵀA{c}) and Frobenius-normalize.
%                   0 → use raw A{c}.
% wINI       : Initial weights on zINI (vector of length numel(zINI)).
% Wfe        : FE integration weights/measures (nIP×1), used for integrals.
% Q          : Columns to build a Wfe-orthonormal projector; expected already
%              orthonormal in the Wfe sense (Qᵀ diag(Wfe) Q = I).
%
% OUTPUTS
% -------
% U          : 1×nC cell; U{c} is the **restricted** basis (size: |zINI|×r_c).
% b          : 1×nC cell; b{c} are target integrals of basis columns (r_c×1).
% wADAPT     : |zINI|×nC matrix; initial adaptive weights (each column = wINI).
% DATAoffline: Returned with defaults filled (idempotent if already set).
% numberMODES: nC×1 vector; r_c = length(b{c}) for each cluster.
%
% ALGORITHM (per cluster c)
% -------------------------
% 1) If Project_SNAP_withQ__MAW_ECM == 1:
%       Aproj = Q(Qᵀ A{c}); scale Aproj to unit Frobenius norm.
%    Else:
%       Aproj = A{c}.
% 2) Orthonormalize columns: Ufull = SVDT([Ufixed, Aproj], TOLloc=0, relative).
% 3) Compute exactness targets against **global** grid:
%       IntExact  = Ufullᵀ * Wfe;
%       IntApprox = Ufull(zINI,:)ᵀ * wINI   % ECM integral with current rule
%    If projecting with Q, set b{c} = IntApprox (consistency with projected
%    subspace); otherwise b{c} = IntExact.
% 4) Restrict basis rows to zINI: U{c} = Ufull(zINI,:).
% 5) Initialize wADAPT(:,c) = wINI and record numberMODES(c) = length(b{c}).
%
% NUMERICAL NOTES
% ---------------
% • Frobenius normalization of the projected snapshots keeps the plastic
%   content at comparable scales before concatenation with Ufixed.
% • Using b{c} = IntApprox when Q-projection is active enforces consistency
%   with the reduced (projected) content actually representable on zINI.
% • SVDT is called with relative scaling (DATAddd.ISRELATIVE=1) and TOLloc=0
%   here; downstream truncation (if any) is handled elsewhere.
%
% ASSUMPTIONS & LIMITATIONS
% -------------------------
% • No cluster overlap is considered (see *_overlap for the alternative).
% • Q is assumed already Wfe-orthonormal; supply Q = SVDT(bsxfun(@times,Qw,1./sqrt(Wfe)))
%   upstream if needed.
% • zINI must index valid global integration points compatible with Wfe size.
%
% COMPLEXITY
% ----------
% Dominated by SVDT on [Ufixed, Aproj] for each cluster; restriction and
% vector products are O(|zINI|·r_c) per cluster.
%
% DEPENDENCIES
% ------------
% DefaultField, SVDT.
% AUTHOR / PLACE / DATE
% ---------------------
% J.A. Hernández Ortega (JAHO) — Barcelona — 21-Sep-2025
% (Comments updated by ChatGPT-5 Thinking)
%==========================================================================


wADAPT = zeros(length(zINI),length(qLATENT)) ;
numberMODES = zeros(length(qLATENT),1) ;
U = cell(size(A));
b = cell(size(A,1),1) ;
DATAoffline = DefaultField(DATAoffline,'Project_SNAP_withQ__MAW_ECM',1) ; % = 1; 

PROJECT_Qapprox  =DATAoffline.Project_SNAP_withQ__MAW_ECM;
DATAddd.ISRELATIVE = 1;
TOLloc =0; %DATAoffline.errorFINT;
for icluster = 1:length(A)
    
    if PROJECT_Qapprox == 1
        Aproy =  Q*(Q'*A{icluster}) ;
        Aproy = Aproy/norm(Aproy,'fro') ;
        %Aproy = SVDT(Aproy)  ;
        
        [ U{icluster},SSS,VVV]  = SVDT([Ufixed,Aproy],TOLloc,DATAddd) ;
        
        
    else
        
        [U{icluster},SSS,VVV]  = SVDT([Ufixed,A{icluster}]) ;
    end
    
    IntExact = U{icluster}'*Wfe ;
    IntApprox = U{icluster}(zINI,:)'*wINI ;
 %   errorINT = norm(IntExact-IntApprox)/norm(IntExact)*100 ;  % Just to check that
    % when PROJECT_Qapprox, errorINT \approx 0
    
    if PROJECT_Qapprox == 1
        b{icluster}  = IntApprox ;   % This is the EXACT integral of this vector
    else
        b{icluster}  = IntExact ;
    end
    
    U{icluster}  = U{icluster}(zINI,:) ;
    %   b{icluster}  = U{icluster}'*wINI ; % This is the ECM integral of this vector
    
    
    wADAPT(:,icluster) = wINI ;
    numberMODES(icluster) = length(b{icluster}) ;
end