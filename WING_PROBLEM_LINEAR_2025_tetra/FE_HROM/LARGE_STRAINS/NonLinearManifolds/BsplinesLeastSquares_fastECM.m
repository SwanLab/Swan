function [DATA_regress_eta_der, nREDcoor] = BsplinesLeastSquares_fastECM(DATA_interp, qINF, VrightT, UU, SS)
%--------------------------------------------------------------------------
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
%     - Produces a compact evaluation structure (`DATA_regress_eta_der`) that
%       allows efficient, vectorized evaluation of η(q), η′(q), η″(q) without
%       creating per-mode anonymous functions.
%
% IMPROVEMENTS OVER PREVIOUS VERSIONS:
%   - Vectorized global spline fit (`spGLO`) for all modes at once.
%   - Precomputes derivatives (`spGLOder1`, `spGLOder2`) and extrapolation
%     boundary data in a single pass.
%   - Stores all data in `DATA_regress_eta_der` for reuse in the online stage,
%     eliminating the performance penalty of repeated spline setup.
%   - Keeps per-mode splines (`sp`, `sp1`, `sp2`) only for plotting and
%     diagnostics.
%
% METHOD OVERVIEW (same theoretical framework as `BsplinesLeastSquares_ECM`):
%   1. Fit cubic B-splines to each mode in `VrightT` with `spap2` (least-squares).
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
%   - qINF        : [1 × N] vector of master-point coordinates
%   - VrightT     : [m × N] matrix of function values (slave modes)
%   - UU          : Projection matrix (e.g., left singular vectors)
%   - SS          : Scaling vector (e.g., singular values)
%
% OUTPUTS:
%   - DATA_regress_eta_der : Struct for downstream evaluation in adaptive ECM:
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

%--------------------------------------------------------------------------



if nargin == 0
    load('tmp1.mat');  % Load fallback data for demo/testing
end

%DATA_interp = DefaultField(DATA_interp, 'INCLUDE_SECOND_DERIVATIVES', 0);
DATA_interp = DefaultField(DATA_interp, 'order_Bslines', 4);
DATA_interp = DefaultField(DATA_interp, 'PortionExtrapolation_plot', 0.2);
L_extrap = DATA_interp.PortionExtrapolation_plot ;
nREDcoor = 1;

DATA_interp = DefaultField(DATA_interp,'Legend_var_x','qINF') ;
DATA_interp = DefaultField(DATA_interp,'Legend_var_y','VrightT') ;
DATA_interp = DefaultField(DATA_interp,'TITLE_additional','') ;

title_ADD = DATA_interp.TITLE_additional ; 

DATA_interp = DefaultField(DATA_interp,'Nfigures_base',15) ;


% Validate inputs
% assert(isvector(qINF), 'qINF must be a vector');
[m, N] = size(VrightT);
% assert(length(qINF) == N, 'qINF and VrightT must have compatible dimensions');
x_min = min(qINF);
x_max = max(qINF) ;
%

% REGRESSION USING qINF as input and VrightT(i,:) as output
% -------------------------------------------------------
%x = qINF ;
order_Bslines = DATA_interp.order_Bslines;
%

knots_number = DATA_interp.NSAMPLES ;
%CoeffREGR_fun = cell()

% Although it can be done all at once, we prefer to make a loop and
% evaluate one by one each coefficient
% CoeffREGR_fun = spap2(knots_number, order_Bslines, qINF, VrightT); % Function
% CoeffREGR_der1 =  fnder(CoeffREGR_fun) ;  % First derivative
% CoeffREGR_der2 =  fnder(CoeffREGR_der1) ;  % Second derivative

% Now let us plot the function and its derivative
sp = cell(size(VrightT,1),1) ;
sp1 = cell(size(VrightT,1),1) ;
sp2 =  cell(size(VrightT,1),1) ;
L = x_max - x_min;

x_plot = linspace(x_min - L_extrap * L, x_max + L_extrap* L, 500);

figure(200)
hold on
xlabel('Amst')
ylabel('Right sing. vect, slave ECM points')
title(['Right sing. vect, slave ECM points + approx. Bsplines',' ',title_ADD])
m = size(VrightT,1) ;
colors = lines(m);
xline(x_min, '--k', 'LineWidth', 1.2, 'DisplayName', 'Training min');
xline(x_max, '--k', 'LineWidth', 1.2, 'DisplayName', 'Training max');

SHOW__derivative = 1;

if  SHOW__derivative == 1
    figure(201  )
    hold on
    xlabel('Amst')
    ylabel('Right sing. vect (1st derivative)')
    title(['ECM slave points: 1st derivative of  approx. Bsplines',' ',title_ADD])
    
    xline(x_min, '--k', 'LineWidth', 1.2, 'DisplayName', 'Training min');
    xline(x_max, '--k', 'LineWidth', 1.2, 'DisplayName', 'Training max');
    
    figure(202)
    hold on
    xlabel('Amst')
    ylabel('Right sing. vect (2nd derivative)')
    title(['ECM slave points:  2nd derivative of approx. Bsplines',' ',title_ADD])
    
    xline(x_min, '--k', 'LineWidth', 1.2, 'DisplayName', 'Training min');
    xline(x_max, '--k', 'LineWidth', 1.2, 'DisplayName', 'Training max');
    
end


% This is for doing the regression for all slaves functions at once
spGLO = spap2(knots_number, order_Bslines, qINF, VrightT);
spGLOder1 =fnder(spGLO) ;
spGLOder2 =fnder(spGLOder1) ;


for islave = 1:size(VrightT,1)
    sp{islave} = spap2(knots_number, order_Bslines, qINF, VrightT(islave,:)); % Function
    sp1{islave} =  fnder(sp{islave}) ;  % First derivative
    sp2{islave} =  fnder(sp1{islave}) ;  % Second derivative
    
    
    [OUTPUT ] = evaluate_spline_with_extrapolation( sp{islave} , sp1{islave} ,  sp2{islave}, x_plot) ;
    figure(200)
    hold on
    plot(qINF, VrightT(islave,:),'DisplayName',['V',num2str(islave)],'Color',colors(islave,:)) ;
    
    plot(x_plot,OUTPUT.VALUE,'--','DisplayName',['Bspline ',num2str(islave)],'Color',colors(islave,:)) ;
    
    if  SHOW__derivative == 1
        figure(201)
        hold on
        plot(x_plot,OUTPUT.DER1,'--','DisplayName',['Bspline ',num2str(islave)],'Color',colors(islave,:)) ;
        figure(202)
        hold on
        plot(x_plot,OUTPUT.DER2,'--','DisplayName',['Bspline ',num2str(islave)],'Color',colors(islave,:)) ;
    end
    
    
end
figure(200)
legend show

figure(201)
legend show

figure(202)
legend show


% USE_VECTORIZED_VERSION =1;
% 
% if USE_VECTORIZED_VERSION == 1
    DATA_regress_eta_der.nameFunctionEvaluate = 'etaFUN_1paramVECT'; % Name of the function
    % INPUTS FOR THE ABOVE FUNCTION
    DATA_regress_eta_der.sp = spGLO ;  % Function
    DATA_regress_eta_der.sp1 = spGLOder1 ; % Derivatives function
    DATA_regress_eta_der.sp2 = spGLOder2 ; % Derivatives function
    DATA_regress_eta_der.xmin = x_min ;
    DATA_regress_eta_der.xmax = x_max ;
    DATA_regress_eta_der.UleftSingular = UU ;
    DATA_regress_eta_der.SSingular = SS ;
    
    % Infor for extrapolation
    INFO_EXTRAP = [] ;
    INFO_EXTRAP.XMIN.s0 = fnval(spGLO,x_min) ;
    INFO_EXTRAP.XMIN.s1 = fnval(spGLOder1,x_min) ;
    INFO_EXTRAP.XMIN.s2 = fnval(spGLOder2,x_min) ;
    INFO_EXTRAP.XMAX.s0 = fnval(spGLO,x_max) ;
    INFO_EXTRAP.XMAX.s1 = fnval(spGLOder1,x_max) ;
    INFO_EXTRAP.XMAX.s2 = fnval(spGLOder2,x_max) ;
    DATA_regress_eta_der.INFO_EXTRAP = INFO_EXTRAP;
    
    disp('Testing etaFUN_1paramVECT')
    qMASTER = [-20,20] ;
    [eta,etaDER1,etaDER2] = etaFUN_1paramVECT(qMASTER,DATA_regress_eta_der)  ;
    
%     
% else
%     DATA_regress_eta_der.nameFunctionEvaluate = 'tauFUN_1param'; % Name of the function
%     % INPUTS FOR THE ABOVE FUNCTION
%     DATA_regress_eta_der.sp = sp ;  % Function
%     DATA_regress_eta_der.sp1 = sp1 ; % Derivatives function
%     DATA_regress_eta_der.sp2 = sp2 ; % Derivatives function
%     DATA_regress_eta_der.xmin = x_min ;
%     DATA_regress_eta_der.xmax = x_max ;
%     DATA_regress_eta_der.UleftSingular = UU ;
%     DATA_regress_eta_der.SSingular = SS ;
%     
%     disp('Testing tauFUN_1param')
%     qMASTER =  [-20,20] ;
%     [tau,tauDER1,tauDER2] = tauFUN_1param(qMASTER,DATA_regress_eta_der)  ;
%     
%     
% end
% 























%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
% BSplines = struct('spline', {}, 'derivative', {}, ...
%     'evaluate', {}, 'evaluate_derivative', {}, 'evaluate_second_derivative', {});
%
% % Setup plots for visualization
% figure(DATA_interp.Nfigures_base); hold on; colors = lines(m);
% xlabel(DATA_interp.Legend_var_x); ylabel(DATA_interp.Legend_var_y); title(['B-spline Approximation of ',DATA_interp.Legend_var_x,...
%     ' vs ', DATA_interp.Legend_var_y]);
%
% xline(x_min, '--k', 'LineWidth', 1.2, 'DisplayName', 'Training min');
% xline(x_max, '--k', 'LineWidth', 1.2, 'DisplayName', 'Training max');
% figure(DATA_interp.Nfigures_base+1); hold on; xlabel('Amst'); ylabel(['First deriv.']);
% title(['B-spline Approximation of 1st deriv',DATA_interp.Legend_var_x,...
%     ' vs ', DATA_interp.Legend_var_y]);
% xline(x_min, '--k', 'LineWidth', 1.2, 'DisplayName', 'Training min');
% xline(x_max, '--k', 'LineWidth', 1.2, 'DisplayName', 'Training max');
% figure(DATA_interp.Nfigures_base+2); hold on; xlabel('Amst'); ylabel('Second deriv.');
% title(['B-spline Approximation of second deriv. of ',DATA_interp.Legend_var_x,...
%     ' vs ', DATA_interp.Legend_var_y]);
% xline(x_min, '--k', 'LineWidth', 1.2, 'DisplayName', 'Training min');
% xline(x_max, '--k', 'LineWidth', 1.2, 'DisplayName', 'Training max');
%
% % Loop over each output function to build its spline approximation
% order_Bslines = DATA_interp.order_Bslines;
% for i = 1:m
%     y = VrightT(i, :);  % i-th function sampled at qINF
%
%     % Least-squares B-spline fit using spap2
%     sp = spap2(DATA_interp.NSAMPLES, order_Bslines, qINF, y);
%     sp1 = fnder(sp);
%     sp2 = fnder(sp1);
%     %  BSplines(i).evaluate = @(x) evaluate_with_quadratic_extrap(sp, sp1, sp2, qINF);
%
%
%     OUTPUT = @(x) evaluate_spline_with_extrapolation(sp, sp1, sp2, x);
%
%
%
%     BSplines(i).evaluate = @(x) OUTPUT(x).VALUE;         % just value
%     BSplines(i).evaluate_derivative  = @(x) OUTPUT(x).DER1;  % will extract dy later
%     BSplines(i).evaluate_second_derivative = @(x) OUTPUT(x).DER2;  % d2y later
%
%
%     %
%     %     sp_der = fnder(sp);
%     %     if DATA_interp.INCLUDE_SECOND_DERIVATIVES == 1
%     %         sp_der_der = fnder(sp_der);
%     %     else
%     %         sp_der_der = [];
%     %     end
%
%     % Store spline and its evaluators
%     %     BSplines(i).spline = sp;
%     %     BSplines(i).derivative = sp_der;
%     %     BSplines(i).evaluate = @(x) fnval(sp, x);
%     %     BSplines(i).evaluate_derivative = @(x) fnval(sp_der, x);
%     %     if ~isempty(sp_der_der)
%     %         BSplines(i).evaluate_second_derivative = @(x) fnval(sp_der_der, x);
%     %     end
%
%     L = x_max - x_min;
%
%     x_plot = linspace(x_min - L_extrap * L, x_max + L_extrap* L, 500);
%
%     % Plot original and fitted values
%     figure(DATA_interp.Nfigures_base);
%     yy = BSplines(i).evaluate(x_plot);
%     plot(qINF, VrightT(i, :), 'Color', colors(i, :), ...
%         'DisplayName', ['\lambda_{', num2str(i), '}']);
%     plot(x_plot, yy, '--', 'Color', colors(i, :), ...
%         'DisplayName', ['B-spline Fit ', num2str(i)]);
%
%     % Plot first derivative
%     figure(DATA_interp.Nfigures_base+1);
%     dydx = BSplines(i).evaluate_derivative(x_plot);
%     plot(x_plot, dydx, '-', 'DisplayName', ['d\lambda_{', num2str(i), '}/dq']);
%
%
%     % Plot second derivative if computed
%     % if ~isempty(sp_der_der)
%     figure(DATA_interp.Nfigures_base+2);
%     d2ydx2 = BSplines(i).evaluate_second_derivative(x_plot);
%     plot(x_plot, d2ydx2, '-', 'DisplayName', ['d^2\lambda_{', num2str(i), '}/dq^2']);
%
% end
%
% % Show legends
% figure(DATA_interp.Nfigures_base); legend show;
% figure(DATA_interp.Nfigures_base+1); legend show;
% figure(DATA_interp.Nfigures_base+2); legend show;
%
% % Create nonlinear evaluation functions
% [NonLinear_q, NonLinear_q_der, NonLinear_q_der2] = make_BSP_evaluator(BSplines, UU, SS);
% end
%
% % -------------------------------------------------------------------------
% function [NonLinear_q, NonLinear_q_der, NonLinear_q_der2] = make_BSP_evaluator(BSplines, UU, SS)
% % make_BSP_evaluator - Creates nonlinear evaluation function handles
% %
% % Inputs:
% %   BSplines - struct array of spline functions and derivatives
% %   UU       - projection matrix (e.g., reduced basis)
% %   SS       - scaling vector for each interpolated function
% %
% % Outputs:
% %   NonLinear_q      - evaluates [x;     UU * (SS .* interp(x))      ]
% %   NonLinear_q_der  - evaluates [1;     UU * (SS .* interp'(x))     ]
% %   NonLinear_q_der2 - evaluates [0;     UU * (SS .* interp''(x))    ]
%
% NonLinear_q      = @(x) fast_evalC(BSplines, x, 'evaluate', UU, SS);
% NonLinear_q_der  = @(x) fast_evalC(BSplines, x, 'evaluate_derivative', UU, SS);
% NonLinear_q_der2 = @(x) fast_evalC(BSplines, x, 'evaluate_second_derivative', UU, SS);
% end
%
% % -------------------------------------------------------------------------
% function OutputGenNON = fast_evalC(BSplines, x, fieldnameLABEL, UU, SS)
% % fast_evalC - Evaluates the spline (or its derivatives) at points x,
% % scales with SS, and projects using UU.
% %
% % Inputs:
% %   BSplines        - struct array with spline evaluators
% %   x               - input values (scalar or vector)
% %   fieldnameLABEL  - 'evaluate', 'evaluate_derivative', or 'evaluate_second_derivative'
% %   UU              - projection matrix (size [nDOF x m])
% %   SS              - scaling vector (1xm)
% %
% % Output:
% %   OutputGenNON - Evaluated and projected result:
% %     - [x; UU * scaled values]               for 'evaluate'
% %     - [1; UU * scaled first derivatives]    for 'evaluate_derivative'
% %     - [0; UU * scaled second derivatives]   for 'evaluate_second_derivative'
%
% x = x(:).';  % Ensure row vector for broadcasting
% m = numel(BSplines);
% n = numel(x);
% coeff = zeros(m, n);
%
% for i = 1:m
%     coeff(i, :) = SS(i) * BSplines(i).(fieldnameLABEL)(x);
% end
%
% switch fieldnameLABEL
%     case 'evaluate'
%         OutputGenNON = [x; UU * coeff];
%     case 'evaluate_derivative'
%         OutputGenNON = [1; UU * coeff];
%     case 'evaluate_second_derivative'
%         OutputGenNON = [0; UU * coeff];
% end
% end
