function [DATA_evaluateTAU_and_DER, nREDcoor] = BsplinesLeastSquares_PLAST_al(DATA_interp, qMASTER_nonl, VrightT, UU, SS)
%--------------------------------------------------------------------------
% This is a modification of BsplinesLeastSquares_PLAST, described below
% The modification is related to the fact that now 
% dL = Phi*tau(q) =  PhiMaster_elast * qMASTER_elast
%            +PhiMaster_nonl*g(qMASTER_plast) ,

%
% BsplinesLeastSquares_PLAST
%
% PURPOSE
%   Fit cubic B-splines (and their derivatives) to represent the nonlinear
%   slave coordinates in a 2-coordinate (elastic–plastic) manifold ROM.
%
%   In this plastic extension, the solution manifold is parameterized by:
%       qL = [ qMASTER_elast , qMASTER_plast ]
%   The displacement decoder has the form:
%       d_L = PhiMaster_elast * qMASTER_elast
%            + PhiMaster_nonl * qMASTER_plast
%            + PhiSlave_nonl  * g(qMASTER_plast) ,
%   where g(·) maps the 1D plastic master coordinate to the reduced slave
%   coordinates. This routine learns g(·) by fitting cubic B-splines to
%   the right singular vectors of the slave coordinate matrix.
%
% INPUTS
%   DATA_interp : struct controlling spline fitting and plotting.
%                 Typical fields:
%       .NSAMPLES                   - number of spline intervals (breaks)
%       .INCLUDE_SECOND_DERIVATIVES - (0 or 1) compute second derivatives
%       .order_Bslines              - spline order (default 4)
%       .PortionExtrapolation_plot  - fraction of domain length to extend in plots
%       .Legend_var_x               - x-axis label for plots
%       .Legend_var_y               - y-axis label for plots
%       .Nfigures_base              - base index for plotting figures
%
%   qMASTER_nonl : [1 × N] vector of plastic master coordinates (training points)
%   VrightT      : [m × N] matrix of function values to fit
%                   (rows = modes in compact SVD coordinates; columns = samples)
%   UU           : [nDOF_slave × m] left singular vectors from compact SVD of slave data
%   SS           : [m × 1] singular values from compact SVD of slave data
%
% OUTPUTS
%   DATA_evaluateTAU_and_DER : struct with all data needed for fast evaluation:
%       .nameFunctionEvaluate - 'tauFUN_1paramELASTplast'
%       .sp, .sp1, .sp2       - vectorized spline (B-form) and its 1st/2nd derivatives
%       .xmin, .xmax          - bounds of qMASTER_plast from training data
%       .UleftSingular        - UU from input
%       .SSingular            - SS from input
%       .INFO_EXTRAP          - boundary values/derivatives for extrapolation
%                               (.XMIN.s0/s1/s2, .XMAX.s0/s1/s2)
%
%   nREDcoor    : number of reduced coordinates defining the manifold (2 for elastic–plastic)
%
% METHOD
%   1) Define the fitting domain [xmin, xmax] from qMASTER_nonl.
%   2) Create global vectorized cubic B-spline fits (`spGLO`) for all rows of VrightT.
%   3) Compute their first and second derivatives (`spGLOder1`, `spGLOder2`).
%   4) For diagnostics, also fit and store per-mode splines (sp{i}, sp1{i}, sp2{i}).
%   5) Plot training data vs. spline fits, plus derivative curves (optional).
%   6) Precompute extrapolation data (function + derivatives at xmin/xmax).
%   7) Package everything into DATA_evaluateTAU_and_DER for use in
%      `tauFUN_1paramELASTplast`, which reconstructs slave contributions in
%      physical space from a given qMASTER_plast.
%
% DEPENDENCIES
%   • spap2   (Curve Fitting Toolbox) for B-spline regression
%   • fnder   for derivative spline construction
%   • fnval   for evaluating B-splines at specific points
%   • evaluate_spline_with_extrapolation (for plotting diagnostics)
%   • tauFUN_1paramELASTplast (for ROM evaluation)
%   • DefaultField (to set default struct fields)
%
% ASSUMPTIONS
%   • qMASTER_nonl is sorted or at least free of duplicates (monotonic for fitting).
%   • VrightT rows correspond to orthonormal SVD modes of slave coordinates.
%   • The elastic coordinate qMASTER_elast is handled separately; this routine
%     only learns g(qMASTER_plast) for the plastic part.
%
% EXAMPLE
%   DATA_interp.NSAMPLES = 20;
%   DATA_interp.INCLUDE_SECOND_DERIVATIVES = 1;
%   [DATA_eval, nRED] = BsplinesLeastSquares_PLAST(DATA_interp, qMASTER_plast, VV', UU, SS);
%   qL = [q_elast; q_plast];
%   [tau, dtau, d2tau] = tauFUN_1paramELASTplast(qL, DATA_eval);
%
% VERSION HISTORY
%   • 2025-07-27 : Adapted from BsplinesLeastSquares_fast for elastic–plastic 2-coordinate manifolds.
%   * Last version: August 13th 2025, Molinos Marfagones, Cartagena.  
%--------------------------------------------------------------------------




if nargin == 0
    load('tmp1.mat');  % Load fallback data for demo/testing
end

DATA_interp = DefaultField(DATA_interp, 'order_Bslines', 4);
DATA_interp = DefaultField(DATA_interp, 'PortionExtrapolation_plot', 0.2);
L_extrap = DATA_interp.PortionExtrapolation_plot ;
nREDcoor = 2;



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

x_plot = linspace(x_min - L_extrap * L, x_max + L_extrap* L, 500);

figure(100  )
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
    figure(101  )
    hold on
    xlabel('qMASTER nonl')
    ylabel('Right sing. vect  (nonlinear slaves) (1st derivative)')
    title('1st derivative of  approx. Bsplines')
    
    xline(x_min, '--k', 'LineWidth', 1.2, 'DisplayName', 'Training min');
    xline(x_max, '--k', 'LineWidth', 1.2, 'DisplayName', 'Training max');
    
    figure(102)
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
    figure(100)
    hold on
    plot(qMASTER_nonl, VrightT(islave,:),'DisplayName',['V',num2str(islave),' S=',num2str(SS(islave)/SS(1))],'Color',colors(islave,:)) ;
    
    plot(x_plot,OUTPUT.VALUE,'--','DisplayName',['Bspline ',num2str(islave)],'Color',colors(islave,:)) ;
    
    if  SHOW__derivative == 1
        figure(101)
        hold on
        plot(x_plot,OUTPUT.DER1,'--','DisplayName',['Bspline ',num2str(islave)],'Color',colors(islave,:)) ;
        figure(102)
        hold on
        plot(x_plot,OUTPUT.DER2,'--','DisplayName',['Bspline ',num2str(islave)],'Color',colors(islave,:)) ;
    end
    
    
end
figure(100)
legend show

figure(101)
legend show

figure(102)
legend show

%
% USE_VECTORIZED_VERSION =1;
%
%     if USE_VECTORIZED_VERSION == 1
DATA_evaluateTAU_and_DER.nameFunctionEvaluate = 'tauFUN_1paramELASTplast'; % Name of the function
% INPUTS FOR THE ABOVE FUNCTION
DATA_evaluateTAU_and_DER.sp = spGLO ;  % Function
DATA_evaluateTAU_and_DER.sp1 = spGLOder1 ; % Derivatives function
DATA_evaluateTAU_and_DER.sp2 = spGLOder2 ; % Derivatives function
DATA_evaluateTAU_and_DER.xmin = x_min ;
DATA_evaluateTAU_and_DER.xmax = x_max ;
DATA_evaluateTAU_and_DER.UleftSingular = UU ;
DATA_evaluateTAU_and_DER.SSingular = SS ;

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
qMASTER = [-0.2,0.2
           -0.2,0.2   ] ;
[tau,tauDER1,tauDER2] = tauFUN_1paramELASTplast(qMASTER,DATA_evaluateTAU_and_DER)  ; 