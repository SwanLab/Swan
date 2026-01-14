function [PhiNONmst,PhiNONslv,qNON,...
    gNONslv,...
    DATA_evaluateTAU_and_DER] = Bsplines_EP_1stmodeDMG(PhiLIN,Phi_non,GmatN,DATA_interp,SNAPdisp_L,GmatN_chol)
% Bsplines_EP_1stmodeDMG
% -------------------------------------------------------------------------
% PURPOSE
%   Construct the *univariate-damage* nonlinear manifold from linear and
%   nonlinear modal snapshots and return:
%     (i)   the *nonlinear master* mode Φ_non,mst and its *slave* complement,
%     (ii)  the master coordinate qNON and the slave amplitudes gNONslv(qNON),
%     (iii) a compact B-spline decoder (and its derivatives) for use in τ(q).
%   This adapts the “first-mode master” pipeline to damage, where the latent
%   parametrization is (qLIN, qNON) and the decoder has the structure
%       τ(q) = [ qLIN ; qLIN*qNON ; qLIN * g(qNON) ].
%
% THEORY (document, pp. 138–146)
%   • One nonlinear master coordinate (qNON) captures the dominant damage
%     evolution; the remaining nonlinear coordinates are treated as *slaves*
%     and regressed as smooth functions of qNON (cubic B-splines).
%   • Casting slave amplitudes as g(qNON) and factoring qLIN yields a
%     differentiable decoder consistent with Newton-type solvers and with
%     hyper-reduction of the manifold equations.
%
% WORKFLOW
%   1) Linear coordinate:
%        qLIN = Φ_linᵀ G SNAP_L.
%   2) Nonlinear raw coordinates:
%        Qnon_all = Φ_nonᵀ G SNAP_L, then normalize by qLIN → Qnon_all./qLIN.
%   3) Monotone master selection:
%        Find a linear combination (ChooseMonotonicLinearCombination) of the
%        normalized nonlinear coordinates that varies *monotonically* with
%        inelastic progression → defines ΦNONmst; normalize in the G-inner
%        product. The corresponding coordinate is qNON.
%   4) Slave basis:
%        Build the G-orthogonal complement of ΦNONmst within Φ_non to obtain
%        ΦNONslv; re-orthogonalize/condition with WSVDT (uses GmatN_chol).
%   5) Slave amplitudes and normalization:
%        qNONslv = ΦNONslvᵀ G SNAP_L, then
%        gNONslv = qNONslv ./ qLIN  (improves regression conditioning).
%        Remove repeated qNON samples to avoid duplicate knots.
%   6) B-spline learning:
%        Fit g(·), g′(·), g″(·) as functions of qNON:
%          – either directly (BsplinesLeastSquares_DAMAGE),
%          – or after uniform subsampling in qNON
%            (BsplinesSUBSAMPL_damage) for robust knot placement.
%   7) Pack evaluators:
%        Return DATA_evaluateTAU_and_DER with spline objects, bounds and
%        extrapolation info so τ(q), ∂τ/∂q, ∂²τ/∂q² can be assembled.
%
% INPUTS
%   PhiLIN       : [nDOF×1] linear/elastic mode (G-orthonormal).
%   Phi_non      : [nDOF×r_nl] nonlinear basis (G-orthonormal columns).
%   GmatN        : [nDOF×nDOF] SPD matrix defining the inner product (e.g., a
%                  constrained stiffness).
%   DATA_interp  : Options for spline fitting and (optional) subsampling:
%                  .NSAMPLES, .order_Bslines, .ratio_NSAMPLES_knots,
%                  .SubSampling_qINELASTIC_master (logical), etc.
%   SNAPdisp_L   : [nDOF×n_snap] snapshot matrix (restricted DOFs).
%   GmatN_chol   : Cholesky-like factor of GmatN for weighted SVD routines.
%
% OUTPUTS
%   PhiNONmst    : [nDOF×1] *nonlinear master* mode.
%   PhiNONslv    : [nDOF×(r_nl−1)] G-orthonormal *nonlinear slave* modes.
%   qNON         : [1×n_eff] master coordinate (after removing duplicates).
%   gNONslv      : [(r_nl−1)×n_eff] slave amplitudes normalized by qLIN,
%                  i.e., gNONslv(qNON).
%   DATA_evaluateTAU_and_DER :
%                  Struct with spline evaluators for g, g′, g″ and metadata
%                  (bounds, extrapolation) to assemble τ and its derivatives.
%
% KEY NOTES
%   • Monotonic master improves identifiability and spline robustness.
%   • All bases and coordinates are computed in the G-inner product.
%   • Normalization by qLIN enforces the factorized decoder form and yields
%     clean first/second derivatives required in consistent tangents.
%   • Optional subsampling along qNON mitigates snapshot clustering and
%     improves numerical conditioning of the spline fit.
%
% DELIVERABLE
%   A compact, differentiable manifold representation ready for
%   tauFUN_1paramDAMAGE and for consistent-tangent assembly in hyper-reduced
%   Newton iterations.
%
% VERSION / AUTHOR
%   23-Oct-2025 — J.A. Hernández (Barcelona). Damage-adapted first-mode
%   master construction with B-spline regression.
% -------------------------------------------------------------------------




if nargin == 0
    load('tmp1.mat')
    %       DATA_interp.NSAMPLES = 300 ;
    %       DATA_interp.ratio_NSAMPLES_knots = 0.8;
    %       DATA_interp.SubSampling_qnonlinear_master = true;
    close all
end

% NumberOfPoints_masterRANGE = size(SNAPdisp_L,2) ;
%
%
% DATA_interp = DefaultField(DATA_interp,'NumberOfPoints_masterRANGE',NumberOfPoints_masterRANGE) ;
%

%  roposed \textbf{``decoder''}
%
% \begin{equation}
% \begin{split}
%  \dL   =     \Phi{}{}  \tauNON(\q) & =  \PhiLIN{}{} \qLIN + \PhiNONmst{}{}\qLIN \qNON +  \PhiNONslv{}{} \qLIN \gNONslv(\qNON)
%  \\& =   \Par{ \PhiMST{}{} \coldos{1}{\qNON}  +  \PhiNONmst{}{} \gNONslv(\qNON)}
%  \\& =  \qLIN  \rowtres{\PhiLIN{}{}}{\PhiNONmst{}{}}{\PhiNONmst{}{} }  \coltres{\qLIN}{\qLIN\qNON}{\qLIN \gNONslv(\qNON)}
%  \\& =   \Phi{}{}   \tauNON(\q)
%  \end{split}
% \end{equation}

% This is the master, nonlinear mode (first mode)
%PhiNONmst =    Phi_non(:,1) ;  % Master nonlinear mode
%PhiNONslv =  Phi_non(:,2:end) ; % Slave modes, nonlinear
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
% Checking that everything is correct (repeat the procedure)
qNON_before_division =  PhiNONmst'*GmatN*SNAPdisp_L ; 
qNON  =  bsxfun(@times,qNON_before_division',1./qLIN')' ;

[qNON,IndicesNonRepeated,bbb ]= unique(qNON) ; % REmove repeated indices


%  
figure(191)
title('Nonlinear latent variable (qNON) as a function of the loading parameter')
hold on
xlabel('Step in inelastic loading')
ylabel('Modal amplitude   nonlinear latent variable')
 
plot(qNON )
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
gNONslv = gNONslv(:,IndicesNonRepeated) ; 

SHOW_before_SVD = 1; 

if SHOW_before_SVD == 1
figure(456)
hold on 
title(['Slave modes amplituted divided by qLIN (gNONslv) as a function of qNON '])
for iii  =1:size(gNONslv,1)
    plot(qNON,gNONslv(iii,:),'DisplayName',['gNONslv_{',num2str(iii),'}'])
end
legend show
end

  

DATA_interp = DefaultField(DATA_interp,'SubSampling_qINELASTIC_master',true) ;

if ~DATA_interp.SubSampling_qINELASTIC_master
    
    [UU,SS,VV] = SVDT(gNONslv) ;
    [DATA_evaluateTAU_and_DER] = BsplinesLeastSquares_DAMAGE(DATA_interp, qNON, VV', UU, SS) ;
    
else
    [DATA_evaluateTAU_and_DER] =  BsplinesSUBSAMPL_damage(qNON,DATA_interp,gNONslv) ;
    
    
    
end


