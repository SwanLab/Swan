function [DATA_evaluateTAU_and_DER,DATA_interp] =  BsplinesSUBSAMPL(qMASTER_nonl,DATA_interp,qSLAVE_nonl)
% =========================================================================
% BSPLINESSUBSAMPL — Uniform master-abscissa subsampling + spline decoder
% =========================================================================
% PURPOSE
%   Build a stable spline decoder τ(q), τ′(q), τ″(q) for plastic SLAVE
%   coordinates by (i) subsampling the MASTER coordinate range to achieve
%   near-uniform coverage, (ii) compressing SLAVE data with a compact SVD,
%   and (iii) fitting least-squares B-splines to the right singular vectors.
%
% WHAT THIS FUNCTION DOES
%   1) Master-axis subsampling:
%        • Create a uniform grid qMASTER_range over [min(qMASTER), max(qMASTER)].
%        • Select nearest snapshots via knnsearch → IndSelectUNI.
%        • Form subsampled sets (qMASTER_nonl_SAMP, qSLAVE_nonl_SAMP).
%   2) Visual diagnostics:
%        • Figure 35: original vs spline-evaluated amplitudes and their first
%          derivatives (helps spot over/under-smoothing).
%   3) Knot budgeting:
%        • knots_max = (#subsamples − order_Bsplines + 1) * ratio_NSAMPLES_knots.
%        • DATA_interp.NSAMPLES is set to knots_max for the downstream fit
%          (so “samples” here effectively means “knots to allocate”).
%   4) Compact SVD on subsampled SLAVE coordinates:
%        • [U,S,V] = SVDT(qSLAVE_nonl_SAMP).
%   5) B-spline regression:
%        • Fit rows of Vᵀ as functions of qMASTER_nonl_SAMP using spap2/fnder
%          (internally via BsplinesLeastSquares_PLAST).
%        • Return callable handles {τ, τ′, τ″} packaged in DATA_evaluateTAU_and_DER.
%
% INPUTS
%   qMASTER_nonl   : [1 × nsnap] MASTER plastic coordinate per snapshot.
%   DATA_interp    : Struct with interpolation/sampling controls:
%                      · NSAMPLES               — target uniform samples on qMASTER
%                      · order_Bslines          — spline polynomial order (p)
%                      · ratio_NSAMPLES_knots   — 0–1 scaling of available knots
%   qSLAVE_nonl    : [r_sl × nsnap] SLAVE plastic coordinates per snapshot.
%
% OUTPUTS
%   DATA_evaluateTAU_and_DER : Struct with function handles to evaluate
%                              τ(q), τ′(q), τ″(q) on demand.
%   DATA_interp              : Possibly updated (sets NSAMPLES = knots_max).
%
% KEY DETAILS / GOTCHAS
%   • Subsampling combats clustering of training points (e.g., near yield),
%     improving spline conditioning and generalization.
%   • The effective number of knots is limited by the number of unique
%     subsampled indices; ensure length(IndSelectUNI) ≥ order_Bslines.
%   • Visual overlays (Figure 35) compare raw amplitudes vs spline evaluation
%     and their derivatives over the full snapshot index domain.
%
% DEPENDENCIES
%   knnsearch (Statistics & ML Toolbox), spap2/fnder (Curve Fitting Toolbox),
%   evaluate_spline_with_extrapolation, SVDT, BsplinesLeastSquares_PLAST,
%   DefaultField.
%
% VERSION / AUTHORSHIP
%   • 18-AUG-2025 — Initial subsampling + spline pipeline. Cartagena.
%   • 07-NOV-2025 — Comments clarified; knot budgeting made explicit; plots
%                   documented. Barcelona.
%   Author: Joaquín Alberto Hernández Ortega (JAHO) — jhortega@cimne.upc.edu
%   (This header was generated automatically by ChatGPT on 07-NOV-2025.)
% =========================================================================

  % Sampling uniform points from the image of  qMASTER_nonl
    qMASTER_range = linspace(min(qMASTER_nonl),max(qMASTER_nonl),DATA_interp.NSAMPLES) ;
    [IndSelect,~] = knnsearch(qMASTER_nonl',qMASTER_range') ;
    IndSelectUNI = unique(IndSelect) ;
    
    qMASTER_nonl_SAMP = qMASTER_nonl(IndSelectUNI) ;
    qSLAVE_nonl_SAMP = qSLAVE_nonl(:,IndSelectUNI) ;
    
    
    figure(35)
    subplot(2,1,1)
    
    hold on
    xlabel('Snapshots')
    ylabel('Modal amplitudes')
    title(['Assessing the effect of subsampling (nsamples = ',num2str(DATA_interp.NSAMPLES),')'])
    
    
    subplot(2,1,2)
    hold on
    xlabel('Snapshots')
    ylabel('Der. Modal amplitudes')
    
    qALL = [qMASTER_nonl;qSLAVE_nonl] ;
    qSUB = [qMASTER_nonl_SAMP;qSLAVE_nonl_SAMP] ;
    colors = lines(size(qSUB,1));
    
    DATA_interp= DefaultField(DATA_interp,'ratio_NSAMPLES_knots',1);
    
    knots_max=  length(IndSelectUNI)-DATA_interp.order_Bslines+1;
    
    knots_max = knots_max*DATA_interp.ratio_NSAMPLES_knots ;
    
    %n_knots = ceil(DATA_interp.ratio_NSAMPLES_knots*length(IndSelectUNI)) ;
    
    for iii  =1:size(qALL,1)
        subplot(2,1,1)
        hold on
        plot( qALL(iii,:),'DisplayName',['qORIG',num2str(iii)],'Color',colors(iii,:)) ;
        
        sp = spap2(knots_max, DATA_interp.order_Bslines, IndSelectUNI, qSUB(iii,:)); % Function
        sp1 =  fnder(sp ) ;  % First derivative
        sp2 =  fnder(sp1) ;  % Second derivative
        
        
        [OUTPUT ] = evaluate_spline_with_extrapolation( sp, sp1 ,  sp2, 1:size(qALL,2)) ;
        %     q_nonlBSPL(islave,:) = OUTPUT.VALUE ;
        %    fff= max(abs(q_nonl(islave,:)));
        %
        %    q_nonlBSPL(islave,:) = OUTPUT.VALUE ;
        %
        
        
        plot(1:size(qALL,2), OUTPUT.VALUE,'--','DisplayName',['qSUB splin',num2str(iii)],'Color',colors(iii,:)) ;
        
        subplot(2,1,2)
        hold on
        plot(1:size(qALL,2), OUTPUT.DER1,'--','DisplayName',['qSUB splin',num2str(iii)],'Color',colors(iii,:)) ;
        
    end
    legend show
    
    % B-splines, it is necessary to sort the entries of qMASTER_nonl is
    % ascending order
    %         [qMASTER_nonl,indORDER] = sort(qMASTER_nonl) ;
    %
    %         figure(126)
    %         hold on
    %         xlabel('Snapshot ordered')
    %         ylabel('Master function (qMASTER nonl)')
    %         plot(qMASTER_nonl)
    %
    %         qSLAVE_nonl = qSLAVE_nonl(:,indORDER)  ;
    %
    % Now we have to establish the mapping qSLAVE_nonl = f(qMASTER_nonl)
    % We do this in an indirect way, by first making
    %% qSLAVE_nonl = UU*diag(SS)*VV^T
    % and then constructing the mapping VV(:,i) = f(qMASTER_nonl)
    %
    %
    
    %
    
    
    [UU,SS,VV] = SVDT(qSLAVE_nonl_SAMP) ;
    
    %   DATA_interp= DefaultField(DATA_interp,'ratio_NSAMPLES_knots',0.9);
    
    DATA_interp.NSAMPLES =knots_max;
    % figure(345)
    % hold on
    % title('Right singular vectors qSLAVE nonl as a function of qMASTER nonl')
    % xlabel('qMASTER nonl')
    % ylabel('qSLAVE nonl')
    %
    %
    % for iii = 1:size(VV,2)
    %     plot(qMASTER_nonl,VV(:,iii),'DisplayName',['\lambda',num2str(iii),'=',num2str(SS(iii)/SS(1))]) ;
    % end
    [DATA_evaluateTAU_and_DER] = BsplinesLeastSquares_PLAST(DATA_interp, qMASTER_nonl_SAMP, VV', UU, SS) ;