function [U,b,wADAPT,DATAoffline,numberMODES] = MAWecmBasisMATRIX_LOCAL_overlap(Ufixed,zINI,qLATENT,A,DATAoffline,wINI,Wfe,Q)
%==========================================================================
% MAWecmBasisMATRIX_LOCAL_overlap
%
% PURPOSE
% -------
% Build per–latent-cluster local exactness bases U{c}, targets b{c}, and
% initial weights wADAPT(:,c)=wINI **with controlled overlap** between
% neighboring clusters. Overlap injects a truncated portion of neighbors’
% snapshots into the local basis to improve continuity/robustness across
% qLATENT while keeping the fixed scaffold Ufixed in all clusters.
%
% WHEN TO USE
% -----------
% Use this instead of *_noverlap when weight functions or internal forces
% should vary smoothly across clusters, or when non-overlapped bases cause
% feasibility/instability in MAW-ECM updates.
%
% INPUTS
% ------
% Ufixed     : Fixed scaffold common to all clusters (e.g., SVDT([Uel, ones])).
% zINI       : Candidate quadrature-point indices (global → local rows).
% qLATENT    : Latent coordinates (only its length is used here for sizing).
% A          : 1×nC cell; A{c} ∈ R^(nIP×m_c), snapshots for cluster c.
% DATAoffline: Options struct. Recognized fields:
%              • NumberOfOverlappingClustersBothDirections (required):
%                    noverlap ≥ 1 → include neighbors [c−noverlap : c+noverlap].
%              • ToleranceOverlappingClusters (default 1e−4):
%                    truncation tol for neighbor content (see SRSVD below).
%              • Project_SNAP_withQ__MAW_ECM (default 1):
%                    1 → project with Q(Qᵀ·) before mixing & SVDT,
%                    0 → use raw snapshots.
% wINI       : Initial weights on zINI (column replicated to all clusters).
% Wfe        : FE integration weights/measures (nIP×1), for integrals.
% Q          : Wfe-orthonormal projector (Qᵀ diag(Wfe) Q = I). If projecting,
%              supply Q already orthonormalized upstream.
%
% OUTPUTS
% -------
% U          : 1×nC cell; U{c} = restricted (|zINI|×r_c) local basis with overlap.
% b          : 1×nC cell; b{c} (r_c×1) target integrals of U{c}’s columns.
% wADAPT     : |zINI|×nC matrix; every column initialized to wINI.
% DATAoffline: Returned with defaults filled (idempotent if already set).
% numberMODES: nC×1; r_c = length(b{c}) for each cluster.
%
% ALGORITHM (for each cluster c)
% ------------------------------
% 1) Neighbor collection:
%       iINI = max(c−noverlap,1),  iFIN = min(c+noverlap,nC),
%       iNEIG = neighbors in [iINI:iFIN] \ {c},  Aneig = cell2mat(A(iNEIG)).
% 2) Optional projection with Q:
%       A_c   = Q(Qᵀ A{c}),   A_nei = Q(Qᵀ Aneig)    (if Project_SNAP_withQ__MAW_ECM==1)
%       else  A_c = A{c},     A_nei = Aneig.
% 3) Overlap by **selective truncation** (SRSVD):
%       [Aover,~,~] = SRSVD({A_c, A_nei}, [0, TOLn], HIDE_OUTPUT=1)
%       • Keeps A_c intact (tol=0) while injecting only the neighbor content
%         above TOLn relative significance → smoother inter-cluster bases.
% 4) Orthonormalize and restrict:
%       Ufull = SVDT([Ufixed, Aover], TOLloc=0, relative scaling);
%       IntExact  = Ufullᵀ Wfe,   IntApprox = Ufull(zINI,:)ᵀ wINI.
%       If projecting with Q → b{c} = IntApprox (consistency with projected
%       content seen on zINI); else b{c} = IntExact.
%       U{c} = Ufull(zINI,:),  wADAPT(:,c) = wINI,  numberMODES(c)=length(b{c}).
%
% NUMERICAL NOTES
% ---------------
% • SRSVD({A_c, A_nei}, [0, TOLn]) preserves the **cluster’s own subspace**
%   exactly while admitting only neighbor directions whose singular values
%   exceed TOLn (relative). This avoids polluting U{c} with noisy/irrelevant
%   neighbor modes yet enforces overlap continuity.
% • With projection on Q, using b{c}=IntApprox ensures exactness targets are
%   consistent with the (projected) content representable by the ECM rule on
%   zINI; without projection, b{c}=IntExact ties to the full FE inner product.
% • DATAddd.ISRELATIVE=1 and TOLloc=0 here; any further truncation is meant
%   to be controlled by TOLn (neighbors) or upstream choices.
%
% ASSUMPTIONS & LIMITATIONS
% -------------------------
% • qLATENT must be ordered so that neighbors [c−noverlap:c+noverlap] are
%   meaningful along the latent path.
% • Q is assumed Wfe-orthonormal already; orthonormalize upstream if needed.
% • Overlap increases r_c; very large noverlap or loose TOLn may inflate the
%   basis and cost downstream MAW-ECM solves.
%
% COMPLEXITY
% ----------
% Per cluster: SVDT on [Ufixed, Aover] dominates. Building Aneig and the
% SRSVD with small noverlap is typically modest compared to SVDT.
%
% DEPENDENCIES
% ------------
% DefaultField, SRSVD, SVDT.
% AUTHOR / PLACE / DATE
% ---------------------
% J.A. Hernández Ortega (JAHO) — Barcelona — 21-Sep-2025
% (Comments updated by ChatGPT-5 Thinking)
%==========================================================================

if nargin == 0
    load('tmp.mat')
end 

DATAoffline = DefaultField(DATAoffline,'ToleranceOverlappingClusters',1e-4) ; % = 1; 
TOLn = DATAoffline.ToleranceOverlappingClusters ; 
noverlap = DATAoffline.NumberOfOverlappingClustersBothDirections ; 
wADAPT = zeros(length(zINI),length(qLATENT)) ;
numberMODES = zeros(length(qLATENT),1) ;
U = cell(size(A));
b = cell(size(A,1),1) ;
DATAoffline = DefaultField(DATAoffline,'Project_SNAP_withQ__MAW_ECM',1) ; % = 1; 

PROJECT_Qapprox  =DATAoffline.Project_SNAP_withQ__MAW_ECM;
DATAddd.ISRELATIVE = 1;
TOLloc =0; %DATAoffline.errorFINT;
DATAnnn.HIDE_OUTPUT = 1; 
for icluster = 1:length(A)
    
     iINI = max(icluster-noverlap,1) ; 
     iFIN = min(icluster+noverlap,length(A)) ; 
     iNEIG = setdiff(iINI:iFIN,icluster); 
     Aneig = cell2mat(A(iNEIG))  ; 
    
    if PROJECT_Qapprox == 1
        Aproy =  Q*(Q'*A{icluster}) ;
        Aneig = Q*(Q'*Aneig) ;
        
        [Aproy,SSSa,VVVa] = SRSVD({Aproy,Aneig},[0,TOLn],DATAnnn) ; 
        %Aproy = Aproy/norm(Aproy,'fro') ;
        %Aproy = SVDT(Aproy)  ;
        
        [ U{icluster},SSS,VVV]  = SVDT([Ufixed,Aproy],TOLloc,DATAddd) ;
        
        
    else
         [Aproy,SSSa,VVVa] = SRSVD({A{icluster},Aneig},[0,TOLn],DATAnnn) ; 
         
         
        [U{icluster},SSS,VVV]  = SVDT([Ufixed,Aproy]) ;
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