function [DATA_regress_w_der, nREDcoor] = BsplinesLeastSquares_MAW_ECM(DATA_interp, q, w)
%--------------------------------------------------------------------------
% This is a modification of BsplinesLeastSquares_fastECM (a function described thoroughly below)
% The goal of the modification is to construct a B-spline approximation of
% ECM weights as a function of the latent variable q
% JAHO, 13-Sept-2025, Saturday, 21:50. Plane Sarajevo-Gerona.
%--------------------------------------------------------------------------
% BsplinesLeastSquares_fastECM
%
% PURPOSE:
%   Improved and specialized version of `BsplinesLeastSquares_fast`,
%   designed for the adaptive Empirical Cubature Method (ECM) workflow.
%
%   It retains the same theoretical basis and methodology as the original
%   `BsplinesLeastSquares_ECM` and `BsplinesLeastSquares_fast` functions,
%   but is streamlined for the master/slave regression step in adaptive ECM.
%
%   In this context, the function:
%     - Receives the right singular vectors (slave modes) from a weighted SVD
%       of ECM slave-point data.
%     - Fits cubic B-splines (and derivatives) to these modes as functions of
%       the master ECM coordinate.
%     - Produces a compact evaluation structure (`DATA_regress_w_der`) that
%       allows efficient, vectorized evaluation of η(q), η′(q), η″(q) without
%       creating per-mode anonymous functions.
%
% IMPROVEMENTS OVER PREVIOUS VERSIONS:
%   - Vectorized global spline fit (`spGLO`) for all modes at once.
%   - Precomputes derivatives (`spGLOder1`, `spGLOder2`) and extrapolation
%     boundary data in a single pass.
%   - Stores all data in `DATA_regress_w_der` for reuse in the online stage,
%     eliminating the performance penalty of repeated spline setup.
%   - Keeps per-mode splines (`sp`, `sp1`, `sp2`) only for plotting and
%     diagnostics.
%
% METHOD OVERVIEW (same theoretical framework as `BsplinesLeastSquares_ECM`):
%   1. Fit cubic B-splines to each mode in `w` with `spap2` (least-squares).
%   2. Compute first and second derivatives via `fnder`.
%   3. Store vectorized fits and derivative splines for fast evaluation.
%   4. Record extrapolation information (value, first derivative, second
%      derivative at both domain ends) for bounded-domain evaluation.
%   5. Optionally produce diagnostic plots for values and derivatives.
%
% INPUTS:
%   - DATA_interp : Struct controlling spline fitting and plotting:
%       • NSAMPLES                   : Number of spline intervals
%       • INCLUDE_SECOND_DERIVATIVES : (0/1) Compute second derivatives
%       • order_Bslines              : Spline order (default = 4)
%       • PortionExtrapolation_plot  : Fraction of domain length to extend in plots
%       • Legend_var_x, Legend_var_y : Labels for plotting
%       • Nfigures_base              : Base figure index for plotting
%   - q        : [1 × N] vector of master-point coordinates
%   - w     : [m × N] matrix of function values (slave modes)
%   - UU          : Projection matrix (e.g., left singular vectors)
%   - SS          : Scaling vector (e.g., singular values)
%
% OUTPUTS:
%   - DATA_regress_w_der : Struct for downstream evaluation in adaptive ECM:
%       • nameFunctionEvaluate : 'etaFUN_1paramVECT'
%       • sp, sp1, sp2         : Vectorized spline and derivatives (B-form)
%       • xmin, xmax           : Domain bounds
%       • UleftSingular, SSingular : Projection and scaling matrices
%       • INFO_EXTRAP          : Boundary value/derivative data
%   - nREDcoor               : Number of reduced coordinates (always 1 here)
%
% REFERENCES (original theory and methodology):
%   - J.A. Hernández (2025), “Digression: on ECM applied to nonlinear manifolds”.
%   - Bravo et al. (2024), SAW-ECM, *International Journal for Numerical Methods
%     in Engineering*, 125(24), e7590.
%
% REMARKS:
%   - The improvements are purely computational; the interpolation and ROM
%     theory are unchanged from the earlier ECM spline-fitting functions.
%   - Intended to be called from `DiscreteECM_adaptWEIGHTSfst` during the
%     offline stage of adaptive ECM.
%
% DATE & PLACE OF MODIFICATION:
%   11-Aug-2025, Molinos Marfagones, Cartagena.
%
% AUTHOR:
%   Joaquín A. Hernández Ortega (UPC/CIMNE)
%   Comments updated by ChatGPT, 12-Aug-2025
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
if nargin == 0
    load('tmp.mat');  % Load fallback data for demo/testing
    close all
    DATA_interp.PortionExtrapolation_plot = 0.2; 
end



 DATA_interp = DefaultField(DATA_interp, 'Extrapolation_Method_WEIGHTS_projection_active_set', 0);


VOL = max(sum(w)); % Volume (or an estimation)
 

%DATA_interp = DefaultField(DATA_interp, 'INCLUDE_SECOND_DERIVATIVES', 0);
DATA_interp = DefaultField(DATA_interp, 'order_Bslines', 4);
DATA_interp = DefaultField(DATA_interp, 'PortionExtrapolation_plot_WEIGHTS', 0.2);
L_extrap = DATA_interp.PortionExtrapolation_plot_WEIGHTS ;
nREDcoor = 1;

DATA_interp = DefaultField(DATA_interp,'Legend_var_x','q') ;
DATA_interp = DefaultField(DATA_interp,'Legend_var_y','Weights') ;
DATA_interp = DefaultField(DATA_interp,'TITLE_additional','') ;

title_ADD = DATA_interp.TITLE_additional ;

DATA_interp = DefaultField(DATA_interp,'Nfigures_base',15) ;


% Validate inputs
% assert(isvector(q), 'q must be a vector');
[m, N] = size(w);
% assert(length(q) == N, 'q and w must have compatible dimensions');
x_min = min(q);
x_max = max(q) ;
%

% REGRESSION USING q as input and w(i,:) as output
% -------------------------------------------------------
%x = q ;
order_Bslines = DATA_interp.order_Bslines;
%

knots_number = DATA_interp.NSAMPLES ;
%CoeffREGR_fun = cell()

% Although it can be done all at once, we prefer to make a loop and
% evaluate one by one each coefficient
% CoeffREGR_fun = spap2(knots_number, order_Bslines, q, w); % Function
% CoeffREGR_der1 =  fnder(CoeffREGR_fun) ;  % First derivative
% CoeffREGR_der2 =  fnder(CoeffREGR_der1) ;  % Second derivative

% Now let us plot the function and its derivative
 
L = x_max - x_min;
NPOINTS_SAMPLE = 1000 ;  
q_plot = linspace(x_min - L_extrap * L, x_max + L_extrap* L, NPOINTS_SAMPLE);

figure(200)
hold on
grid on
xlabel('qLATENT')
ylabel(' weights ECM points div. by total VOLUME')
title([' weights ECM points + approx. Bsplines, div. by total VOLUME',' ',title_ADD])
m = size(w,1) ;
colors = lines(m);
xline(x_min, '--k', 'LineWidth', 1.2, 'DisplayName', 'Training min');
xline(x_max, '--k', 'LineWidth', 1.2, 'DisplayName', 'Training max');

SHOW__derivative = 1;

if  SHOW__derivative == 1
    figure(201  )
    hold on
    xlabel('qLATENT')
    ylabel('1st. der. weights ECM points')
    title(['1st. der  weights ECM points ( approx. Bsplines)',' ',title_ADD])
    
    xline(x_min, '--k', 'LineWidth', 1.2, 'DisplayName', 'Training min');
    xline(x_max, '--k', 'LineWidth', 1.2, 'DisplayName', 'Training max');
    
    figure(202)
    hold on
    xlabel('q')
    ylabel('2nd der.   weights ECM points')
    title(['2nd der weights ECM points ( approx. Bsplines)',' ',title_ADD])
    
    xline(x_min, '--k', 'LineWidth', 1.2, 'DisplayName', 'Training min');
    xline(x_max, '--k', 'LineWidth', 1.2, 'DisplayName', 'Training max');
    
end


% This is for doing the regression for all slaves functions at once


spGLO = spap2(knots_number, order_Bslines, q, w);
spGLOder1 =fnder(spGLO) ;
spGLOder2 =fnder(spGLOder1) ;


% USE_VECTORIZED_VERSION =1;
%
% if USE_VECTORIZED_VERSION == 1


% INPUTS FOR THE ABOVE FUNCTION
DATA_regress_w_der.sp = spGLO ;  % Function
DATA_regress_w_der.sp1 = spGLOder1 ; % Derivatives function
DATA_regress_w_der.sp2 = spGLOder2 ; % Derivatives function
DATA_regress_w_der.xmin = x_min ;
DATA_regress_w_der.xmax = x_max ;
 
% Infor for extrapolation
INFO_EXTRAP = [] ;
INFO_EXTRAP.XMIN.s0 = fnval(spGLO,x_min) ;
INFO_EXTRAP.XMIN.s1 = fnval(spGLOder1,x_min) ;
INFO_EXTRAP.XMIN.s2 = fnval(spGLOder2,x_min) ;
INFO_EXTRAP.XMAX.s0 = fnval(spGLO,x_max) ;
INFO_EXTRAP.XMAX.s1 = fnval(spGLOder1,x_max) ;
INFO_EXTRAP.XMAX.s2 = fnval(spGLOder2,x_max) ;
DATA_regress_w_der.INFO_EXTRAP = INFO_EXTRAP;
DATA_regress_w_der.VOL = VOL ; 

%disp('Testing etaFUN_1paramVECT')
 

if DATA_interp.Extrapolation_Method_WEIGHTS_projection_active_set == 0

DATA_regress_w_der.nameFunctionEvaluate = 'wFUN_1paramVECT'; % Name of the function
[wAPPROX,wAPPROX_DER1,wAPPROX_DER2] = wFUN_1paramVECT(q_plot,DATA_regress_w_der)  ;

else
    DATA_regress_w_der.nameFunctionEvaluate = 'wFUN_1paramVECTp';
    [wAPPROX,wAPPROX_DER1,wAPPROX_DER2] = wFUN_1paramVECTp(q_plot,DATA_regress_w_der)  ;

    
end


 


% DATA_interp.setElements =setElements_cand ;



for ipoint = 1:size(w,1)
     figure(200)
    hold on
    plot(q, w(ipoint,:)/VOL,'DisplayName',['wREL',num2str(ipoint),' e = ',num2str(DATA_interp.setElements(ipoint))],'Color',colors(ipoint,:)) ;
    
    plot(q_plot,wAPPROX(ipoint,:)/VOL,'--','DisplayName',['Bspline ',num2str(ipoint)],'Color',colors(ipoint,:)) ;
    
    if  SHOW__derivative == 1
        figure(201)
        hold on
        plot(q_plot,wAPPROX_DER1(ipoint,:),'--','DisplayName',['Bspline ',num2str(ipoint)],'Color',colors(ipoint,:)) ;
        figure(202)
        hold on
        plot(q_plot,wAPPROX_DER2(ipoint,:),'--','DisplayName',['Bspline ',num2str(ipoint)],'Color',colors(ipoint,:)) ;
    end
    
    
end
figure(200)
legend show

figure(201)
legend show

figure(202)
legend show


 