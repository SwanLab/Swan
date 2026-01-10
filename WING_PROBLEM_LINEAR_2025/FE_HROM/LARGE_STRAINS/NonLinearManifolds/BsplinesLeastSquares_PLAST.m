function [DATA_evaluateTAU_and_DER] = BsplinesLeastSquares_PLAST(DATA_interp, qMASTER_nonl, VrightT, UU, SS)
% =========================================================================
% BSPLINESLEASTSQUARES_PLAST — LSQ B-splines for plastic slave decoder τ(q)
% =========================================================================
% PURPOSE
%   Fit least-squares B-splines (and their derivatives) to represent the
%   NONLINEAR plastic SLAVE coordinates in a 2-coordinate (elastic–plastic)
%   manifold ROM. The manifold is parameterized by qL = [q_ELAST ; q_PLAST],
%   and the displacement decoder reads:
%       d_L(q) = Φ_ELAST*q_ELAST + Φ_master*q_PLAST + Φ_slave * g(q_PLAST),
%   where g(·) is learned here by fitting B-splines to the compact-SVD
%   coordinates (right singular vectors) of the SLAVE data.
%
% WHAT THIS ROUTINE BUILDS
%   • Global, vectorized spline fit for all SLAVE modes at once:
%       spGLO, spGLOder1, spGLOder2  (B-form and its 1st/2nd derivatives)
%   • Per-mode diagnostic splines (optional plots).
%   • Extrapolation metadata at the training interval endpoints
%     (value/1st/2nd derivatives at xmin/xmax) for stable online use.
%   • A callable package DATA_evaluateTAU_and_DER consumed by
%     tauFUN_1paramELASTplast (or tauFUN_1paramELASTplastSYM).
%
% INPUTS
%   DATA_interp : Struct controlling the spline regression and plots
%       .NSAMPLES                  — nominal #knots / intervals for spap2
%       .INCLUDE_SECOND_DERIVATIVES— {0,1} compute/store 2nd derivatives
%       .order_Bslines             — spline order (default 4 = cubic)
%       .PortionExtrapolation_plot — fraction of domain used to extend x-axis
%       .Legend_var_x, .Legend_var_y, .Nfigures_base — plotting cosmetics
%       .ExploitSymmetryPlasticLatentVariable — if true, symmetric τ variant
%
%   qMASTER_nonl : [1 × N] training samples of the MASTER plastic coordinate
%   VrightT      : [m × N] values to fit (right-singular vectors in SLAVE SVD basis)
%   UU           : [n_slave × m] left singular vectors of SLAVE data
%   SS           : [m × 1] singular values of SLAVE data
%
% OUTPUTS
%   DATA_evaluateTAU_and_DER : Struct with all data for fast online eval
%       .nameFunctionEvaluate   — 'tauFUN_1paramELASTplast' or '...SYM'
%       .sp, .sp1, .sp2         — vectorized spline and its 1st/2nd derivs
%       .xmin, .xmax            — training bounds for q_PLAST
%       .UleftSingular, .SSingular — (UU,SS) for compact reconstruction
%       .IndexLatentPlasticVariable — index of q_PLAST in qL
%       .INFO_EXTRAP.XMIN/XMAX.s0/s1/s2 — endpoint data for extrapolation
%
% METHOD (overview)
%   1) Determine training interval: [xmin,xmax] = bounds(qMASTER_nonl).
%   2) Global vectorized fits:
%        spGLO      = spap2(NSAMPLES, order, qMASTER_nonl, VrightT)
%        spGLOder1  = fnder(spGLO);   spGLOder2 = fnder(spGLOder1)
%   3) (Diagnostics) For each row i of VrightT, also store sp{i}, sp1{i}, sp2{i}
%      and plot raw vs spline/derivative curves over an extended domain.
%   4) Package everything into DATA_evaluateTAU_and_DER for online τ(q), τ′(q), τ″(q).
%
% ASSUMPTIONS / TIPS
%   • qMASTER_nonl should be monotone (or at least deduplicated) for stable
%     univariate spline fitting; upstream routines handle sorting/subsampling.
%   • VrightT rows correspond to an orthonormal SVD basis of SLAVE data;
%     using the compact-SVD coordinates improves conditioning and reduces
%     the number of spline regressions.
%   • Choose NSAMPLES ≥ order_Bslines and within the number of unique
%     master samples; increase NSAMPLES to capture sharper plastic transients.
%
% DEPENDENCIES
%   spap2, fnder, fnval (Curve Fitting Toolbox)
%   evaluate_spline_with_extrapolation
%   tauFUN_1paramELASTplast / tauFUN_1paramELASTplastSYM
%   DefaultField
%
% EXAMPLE
%   DATA_interp.NSAMPLES = 20;
%   DATA_interp.INCLUDE_SECOND_DERIVATIVES = 1;
%   [DATA_eval] = BsplinesLeastSquares_PLAST(DATA_interp, qMASTER_plast, VV', UU, SS);
%   [tau, dtau, d2tau] = tauFUN_1paramELASTplast([q_elast; q_plast], DATA_eval);
%
% VERSION HISTORY
%   • 2025-07-27 — Adapted from BsplinesLeastSquares_fast for elastoplastic
%                  2-coordinate manifolds.
%   • 13-AUG-2025 — Last version (previous header). Molinos Marfagones, Cartagena.
%   • 07-NOV-2025 — Comments refreshed; clarified vectorized vs per-mode fits,
%                    endpoint extrapolation, and online packaging. Barcelona.
%   Author: Joaquín Alberto Hernández Ortega (JAHO) — jhortega@cimne.upc.edu
%   (This header was generated automatically by ChatGPT on 07-NOV-2025.)
% =========================================================================
 
%--------------------------------------------------------------------------




if nargin == 0
    load('tmp1.mat');  % Load fallback data for demo/testing
end

DATA_interp = DefaultField(DATA_interp, 'order_Bslines', 4);
DATA_interp = DefaultField(DATA_interp, 'PortionExtrapolation_plot', 0.2);
L_extrap = DATA_interp.PortionExtrapolation_plot ;
%nREDcoor = 2;



% Validate inputs
% assert(isvector(qMASTER_nonl), 'qMASTER_nonl must be a vector');
[m, N] = size(VrightT);
% assert(length(qMASTER_nonl) == N, 'qMASTER_nonl and VrightT must have compatible dimensions');
x_min = min(qMASTER_nonl);
x_max = max(qMASTER_nonl) ;
%

% REGRESSION USING qMASTER_nonl as input and VrightT(i,:) as output
% -------------------------------------------------------
%x = qMASTER_nonl ;
order_Bslines = DATA_interp.order_Bslines;
%

knots_number = DATA_interp.NSAMPLES ;
%CoeffREGR_fun = cell()

% Although it can be done all at once, we prefer to make a loop and
% evaluate one by one each coefficient
% CoeffREGR_fun = spap2(knots_number, order_Bslines, qMASTER_nonl, VrightT); % Function
% CoeffREGR_der1 =  fnder(CoeffREGR_fun) ;  % First derivative
% CoeffREGR_der2 =  fnder(CoeffREGR_der1) ;  % Second derivative

% Now let us plot the function and its derivative
sp = cell(size(VrightT,1),1) ;
sp1 = cell(size(VrightT,1),1) ;
sp2 =  cell(size(VrightT,1),1) ;
L = x_max - x_min;

x_plot = linspace(x_min - L_extrap * L, x_max + L_extrap* L, size(VrightT,2));

figure(100  )
subplot(2,1,1)
hold on
xlabel('qMASTER nonl')
ylabel('Right sing. vect')
title('Right singular vectors (nonlinear slaves) + approx. Bsplines')
m = size(VrightT,1) ;
colors = lines(m);
xline(x_min, '--k', 'LineWidth', 1.2, 'DisplayName', 'Training min');
xline(x_max, '--k', 'LineWidth', 1.2, 'DisplayName', 'Training max');

SHOW__derivative = 1;

if  SHOW__derivative == 1
    figure(100  )
    subplot(2,1,2)

    hold on
    xlabel('qMASTER nonl')
    ylabel('Right sing. vect  (nonlinear slaves) (1st derivative)')
    title('1st derivative of  approx. Bsplines')
    
    xline(x_min, '--k', 'LineWidth', 1.2, 'DisplayName', 'Training min');
    xline(x_max, '--k', 'LineWidth', 1.2, 'DisplayName', 'Training max');
    
    figure(101)
    hold on
    xlabel('qMASTER_nonl')
    ylabel('Right sing. vect (nonlinear slaves) (2nd derivative)')
    title('2nd derivative of approx. Bsplines')
    
    xline(x_min, '--k', 'LineWidth', 1.2, 'DisplayName', 'Training min');
    xline(x_max, '--k', 'LineWidth', 1.2, 'DisplayName', 'Training max');
    
end


% This is for doing the regression for all slaves functions at once
spGLO = spap2(knots_number, order_Bslines, qMASTER_nonl, VrightT);
spGLOder1 =fnder(spGLO) ;
spGLOder2 =fnder(spGLOder1) ;


for islave = 1:size(VrightT,1)
    sp{islave} = spap2(knots_number, order_Bslines, qMASTER_nonl, VrightT(islave,:)); % Function
    sp1{islave} =  fnder(sp{islave}) ;  % First derivative
    sp2{islave} =  fnder(sp1{islave}) ;  % Second derivative
    
    
    [OUTPUT ] = evaluate_spline_with_extrapolation( sp{islave} , sp1{islave} ,  sp2{islave}, x_plot) ;
     figure(100  )
    subplot(2,1,1)

    hold on
    plot(qMASTER_nonl, VrightT(islave,:),'DisplayName',['V',num2str(islave),' S=',num2str(SS(islave)/SS(1))],'Color',colors(islave,:)) ;
    
    plot(x_plot,OUTPUT.VALUE,'--','DisplayName',['Bspline ',num2str(islave)],'Color',colors(islave,:)) ;
    
    if  SHOW__derivative == 1
         figure(100  )
    subplot(2,1,2)

        hold on
        plot(x_plot,OUTPUT.DER1,'DisplayName',['Bspline ',num2str(islave)],'Color',colors(islave,:)) ;
        figure(101)
        hold on
        plot(x_plot,OUTPUT.DER2,'DisplayName',['Bspline ',num2str(islave)],'Color',colors(islave,:)) ;
    end
    
    
end
figure(100)
legend show

figure(101)
legend show
 

%
% USE_VECTORIZED_VERSION =1;
%
%     if USE_VECTORIZED_VERSION == 1
if DATA_interp.ExploitSymmetryPlasticLatentVariable
    DATA_evaluateTAU_and_DER.nameFunctionEvaluate = 'tauFUN_1paramELASTplastSYM';
    DATA_evaluateTAU_and_DER.signTRAIN = DATA_interp.signTRAIN ; 
else
DATA_evaluateTAU_and_DER.nameFunctionEvaluate = 'tauFUN_1paramELASTplast'; % Name of the function
end


% INPUTS FOR THE ABOVE FUNCTION
DATA_evaluateTAU_and_DER.sp = spGLO ;  % Function
DATA_evaluateTAU_and_DER.sp1 = spGLOder1 ; % Derivatives function
DATA_evaluateTAU_and_DER.sp2 = spGLOder2 ; % Derivatives function
DATA_evaluateTAU_and_DER.xmin = x_min ;
DATA_evaluateTAU_and_DER.xmax = x_max ;
DATA_evaluateTAU_and_DER.UleftSingular = UU ;
DATA_evaluateTAU_and_DER.SSingular = SS ;
DATA_interp = DefaultField(DATA_interp,'IndexLatentPlasticVariable',2) ; 
DATA_evaluateTAU_and_DER.IndexLatentPlasticVariable = DATA_interp.IndexLatentPlasticVariable ; 
% Infor for extrapolation
INFO_EXTRAP = [] ;
INFO_EXTRAP.XMIN.s0 = fnval(spGLO,x_min) ;
INFO_EXTRAP.XMIN.s1 = fnval(spGLOder1,x_min) ;
INFO_EXTRAP.XMIN.s2 = fnval(spGLOder2,x_min) ;
INFO_EXTRAP.XMAX.s0 = fnval(spGLO,x_max) ;
INFO_EXTRAP.XMAX.s1 = fnval(spGLOder1,x_max) ;
INFO_EXTRAP.XMAX.s2 = fnval(spGLOder2,x_max) ;
DATA_evaluateTAU_and_DER.INFO_EXTRAP = INFO_EXTRAP;

disp('Testing tauFUN_1paramELASTplast')

if DATA_interp.ExploitSymmetryPlasticLatentVariable
   disp('Testing tauFUN_1paramELASTplastsYM')
   qMASTER = zeros(DATA_interp.nREDcoor,2) ; 
qMASTER(DATA_interp.IndexLatentPlasticVariable,:) = [ 0.2,0.3] ; 
[tau,tauDER1,tauDER2] = tauFUN_1paramELASTplastSYM(qMASTER,DATA_evaluateTAU_and_DER)  ; 



else
disp('Testing tauFUN_1paramELASTplast')
qMASTER = zeros(DATA_interp.nREDcoor,2) ; 
qMASTER(DATA_interp.IndexLatentPlasticVariable,:) = [-0.2,0.2] ; 
[tau,tauDER1,tauDER2] = tauFUN_1paramELASTplast(qMASTER,DATA_evaluateTAU_and_DER)  ; 



end
% qMASTER = [0,0.2
%            -0.2,0.2   ] ;




 