function [NonLinear_q, NonLinear_q_der, NonLinear_q_der2, nREDcoor]= BsplinesLeastSquares_ECMfst(DATA_interp, qINF, qSUP, UU, SS)
% This is a modification of BsplinesLeastSquares_ECM (described below) for
% handling more efficiently the regression tasks ---as done in
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/LARGE_STRAINS/NonLinearManifolds/BsplinesLeastSquares_fast.m
% JAHO, 11th August 2025, Molinos MArfagones, Cartagena.

%--------------------------------------------------------------------------
% function [NonLinear_q, NonLinear_q_der, NonLinear_q_der2, nREDcoor] = ...
%     BsplinesLeastSquares(DATA_interp, qINF, qSUP, UU, SS)
%
% PURPOSE:
%   Constructs nonlinear mappings based on B-spline least-squares fitting
%   of a set of output functions sampled over a 1D input domain.
%   These mappings are useful in reduced-order modeling (ROM) or nonlinear
%   parameterizations where a low-dimensional input variable (qINF) governs
%   a higher-dimensional response (qSUP).
%
% FUNCTIONALITY:
%   - Each row of `qSUP` is approximated by a univariate B-spline function
%     fitted to data points (qINF, qSUP(i,:)) using `spap2`.
%   - Derivatives are computed via `fnder` and can be evaluated up to second order.
%   - Evaluation outside the training domain is handled with **quadratic extrapolation**
%     using a Taylor expansion.
%   - The final outputs are function handles representing:
%       * NonLinear_q(x)      = [x; UU * (SS .* interp(x))]
%       * NonLinear_q_der(x)  = [1; UU * (SS .* interp'(x))]
%       * NonLinear_q_der2(x) = [0; UU * (SS .* interp''(x))]
%
% INPUTS:
%   - DATA_interp : struct with the following fields:
%       * NSAMPLES                    - number of subintervals (B-spline breaks)
%       * INCLUDE_SECOND_DERIVATIVES - flag to compute 2nd derivatives (0 or 1)
%       * order_Bslines              - (optional) spline order (default = 4)
%       * PortionExtrapolation_plot  - range of extrapolation to show in plots
%       * Legend_var_x               - label for x-axis in plots (default: 'qINF')
%       * Legend_var_y               - label for y-axis in plots (default: 'qSUP')
%       * Nfigures_base              - figure index base for plotting
%
%   - qINF : 1xN vector of input coordinates (e.g., reduced variable samples)
%   - qSUP : mxN matrix of output values corresponding to each input (each row is a function)
%   - UU   : projection matrix (size: [nDOF x m]), e.g., reduced basis for outputs
%   - SS   : scaling vector (1 x m), used to scale each interpolated function
%
% OUTPUTS:
%   - NonLinear_q      : handle to nonlinear mapping q ↦ output
%   - NonLinear_q_der  : handle to first derivative of the mapping
%   - NonLinear_q_der2 : handle to second derivative (if enabled)
%   - nREDcoor         : dummy output (always 1), used for legacy compatibility
%
% FEATURES:
%   - Spline fitting performed using `spap2` (least-squares approach).
%   - First and second derivatives obtained using `fnder`.
%   - Quadratic extrapolation used for evaluating points outside the training domain.
%   - Generates diagnostic plots in three figures:
%       * Fig. N        : fitted splines over input domain
%       * Fig. N + 1    : first derivatives of the splines
%       * Fig. N + 2    : second derivatives of the splines (if computed)
%
% USAGE EXAMPLE:
%   DATA_interp.NSAMPLES = 15;
%   DATA_interp.INCLUDE_SECOND_DERIVATIVES = 1;
%   DATA_interp.order_Bslines = 4;  % optional
%   [NLq, NLq_d, NLq_dd] = BsplinesLeastSquares(DATA_interp, qINF, qSUP, UU, SS);
%   y = NLq(0.5);  % Evaluate at q = 0.5
%
% DEPENDENCIES:
%   - Uses `evaluate_spline_with_extrapolation.m` for evaluation logic
%   - Uses `DefaultField.m` to assign defaults
%
% AUTHOR:
%   J.A. Hernández Ortega (UPC/CIMNE), July 2025
%--------------------------------------------------------------------------



if nargin == 0
    load('tmp1.mat');  % Load fallback data for demo/testing
end

DATA_interp = DefaultField(DATA_interp, 'INCLUDE_SECOND_DERIVATIVES', 0);
DATA_interp = DefaultField(DATA_interp, 'order_Bslines', 4);
DATA_interp = DefaultField(DATA_interp, 'PortionExtrapolation_plot', 0.2);
L_extrap = DATA_interp.PortionExtrapolation_plot ;
nREDcoor = 1;

DATA_interp = DefaultField(DATA_interp,'Legend_var_x','qINF') ;
DATA_interp = DefaultField(DATA_interp,'Legend_var_y','qSUP') ;
DATA_interp = DefaultField(DATA_interp,'Nfigures_base',15) ;


% Validate inputs
% assert(isvector(qINF), 'qINF must be a vector');
[m, N] = size(qSUP);
% assert(length(qINF) == N, 'qINF and qSUP must have compatible dimensions');
x_min = min(qINF);
x_max = max(qINF) ;
% Initialize storage for spline approximations and their derivatives
BSplines = struct('spline', {}, 'derivative', {}, ...
    'evaluate', {}, 'evaluate_derivative', {}, 'evaluate_second_derivative', {});

% Setup plots for visualization
figure(DATA_interp.Nfigures_base); hold on; colors = lines(m);
xlabel(DATA_interp.Legend_var_x); ylabel(DATA_interp.Legend_var_y); title(['B-spline Approximation of ',DATA_interp.Legend_var_x,...
    ' vs ', DATA_interp.Legend_var_y]);

xline(x_min, '--k', 'LineWidth', 1.2, 'DisplayName', 'Training min');
xline(x_max, '--k', 'LineWidth', 1.2, 'DisplayName', 'Training max');
figure(DATA_interp.Nfigures_base+1); hold on; xlabel('qINF'); ylabel(['First deriv.']);
title(['B-spline Approximation of 1st deriv',DATA_interp.Legend_var_x,...
    ' vs ', DATA_interp.Legend_var_y]);
xline(x_min, '--k', 'LineWidth', 1.2, 'DisplayName', 'Training min');
xline(x_max, '--k', 'LineWidth', 1.2, 'DisplayName', 'Training max');
figure(DATA_interp.Nfigures_base+2); hold on; xlabel('qINF'); ylabel('Second deriv.');
title(['B-spline Approximation of second deriv. of ',DATA_interp.Legend_var_x,...
    ' vs ', DATA_interp.Legend_var_y]);
xline(x_min, '--k', 'LineWidth', 1.2, 'DisplayName', 'Training min');
xline(x_max, '--k', 'LineWidth', 1.2, 'DisplayName', 'Training max');

% Loop over each output function to build its spline approximation
order_Bslines = DATA_interp.order_Bslines;
for i = 1:m
    y = qSUP(i, :);  % i-th function sampled at qINF
    
    % Least-squares B-spline fit using spap2
    sp = spap2(DATA_interp.NSAMPLES, order_Bslines, qINF, y);
    sp1 = fnder(sp);
    sp2 = fnder(sp1);
    %  BSplines(i).evaluate = @(x) evaluate_with_quadratic_extrap(sp, sp1, sp2, qINF);
    
    
    OUTPUT = @(x) evaluate_spline_with_extrapolation(sp, sp1, sp2, x);
    
    
    
    BSplines(i).evaluate = @(x) OUTPUT(x).VALUE;         % just value
    BSplines(i).evaluate_derivative  = @(x) OUTPUT(x).DER1;  % will extract dy later
    BSplines(i).evaluate_second_derivative = @(x) OUTPUT(x).DER2;  % d2y later
    
    
    %
    %     sp_der = fnder(sp);
    %     if DATA_interp.INCLUDE_SECOND_DERIVATIVES == 1
    %         sp_der_der = fnder(sp_der);
    %     else
    %         sp_der_der = [];
    %     end
    
    % Store spline and its evaluators
    %     BSplines(i).spline = sp;
    %     BSplines(i).derivative = sp_der;
    %     BSplines(i).evaluate = @(x) fnval(sp, x);
    %     BSplines(i).evaluate_derivative = @(x) fnval(sp_der, x);
    %     if ~isempty(sp_der_der)
    %         BSplines(i).evaluate_second_derivative = @(x) fnval(sp_der_der, x);
    %     end
    
    L = x_max - x_min;
    
    x_plot = linspace(x_min - L_extrap * L, x_max + L_extrap* L, 500);
    
    % Plot original and fitted values
    figure(DATA_interp.Nfigures_base);
    yy = BSplines(i).evaluate(x_plot);
    plot(qINF, qSUP(i, :), 'Color', colors(i, :), ...
        'DisplayName', ['\lambda_{', num2str(i), '}']);
    plot(x_plot, yy, '--', 'Color', colors(i, :), ...
        'DisplayName', ['B-spline Fit ', num2str(i)]);
    
    % Plot first derivative
    figure(DATA_interp.Nfigures_base+1);
    dydx = BSplines(i).evaluate_derivative(x_plot);
    plot(x_plot, dydx, '-', 'DisplayName', ['d\lambda_{', num2str(i), '}/dq']);
    
    
    % Plot second derivative if computed
    % if ~isempty(sp_der_der)
    figure(DATA_interp.Nfigures_base+2);
    d2ydx2 = BSplines(i).evaluate_second_derivative(x_plot);
    plot(x_plot, d2ydx2, '-', 'DisplayName', ['d^2\lambda_{', num2str(i), '}/dq^2']);
    
end

% Show legends
figure(DATA_interp.Nfigures_base); legend show;
figure(DATA_interp.Nfigures_base+1); legend show;
figure(DATA_interp.Nfigures_base+2); legend show;

% Create nonlinear evaluation functions
[NonLinear_q, NonLinear_q_der, NonLinear_q_der2] = make_BSP_evaluatorECM(BSplines, UU, SS);
end

% -------------------------------------------------------------------------
function [NonLinear_q, NonLinear_q_der, NonLinear_q_der2] = make_BSP_evaluatorECM(BSplines, UU, SS)
% make_BSP_evaluator - Creates nonlinear evaluation function handles
%
% Inputs:
%   BSplines - struct array of spline functions and derivatives
%   UU       - projection matrix (e.g., reduced basis)
%   SS       - scaling vector for each interpolated function
%
% Outputs:
%   NonLinear_q      - evaluates [x;     UU * (SS .* interp(x))      ]
%   NonLinear_q_der  - evaluates [1;     UU * (SS .* interp'(x))     ]
%   NonLinear_q_der2 - evaluates [0;     UU * (SS .* interp''(x))    ]

NonLinear_q      = @(x) fast_evalC_ECM(BSplines, x, 'evaluate', UU, SS);
NonLinear_q_der  = @(x) fast_evalC_ECM(BSplines, x, 'evaluate_derivative', UU, SS);
NonLinear_q_der2 = @(x) fast_evalC_ECM(BSplines, x, 'evaluate_second_derivative', UU, SS);
end

% -------------------------------------------------------------------------
function OutputGenNON = fast_evalC_ECM(BSplines, x, fieldnameLABEL, UU, SS)
% fast_evalC - Evaluates the spline (or its derivatives) at points x,
% scales with SS, and projects using UU.
%
% Inputs:
%   BSplines        - struct array with spline evaluators
%   x               - input values (scalar or vector)
%   fieldnameLABEL  - 'evaluate', 'evaluate_derivative', or 'evaluate_second_derivative'
%   UU              - projection matrix (size [nDOF x m])
%   SS              - scaling vector (1xm)
%
% Output:
%   OutputGenNON - Evaluated and projected result:
%     - [x; UU * scaled values]               for 'evaluate'
%     - [1; UU * scaled first derivatives]    for 'evaluate_derivative'
%     - [0; UU * scaled second derivatives]   for 'evaluate_second_derivative'

x = x(:).';  % Ensure row vector for broadcasting
m = numel(BSplines);
n = numel(x);
coeff = zeros(m, n);

for i = 1:m
    coeff(i, :) = SS(i) * BSplines(i).(fieldnameLABEL)(x);
end

% switch fieldnameLABEL
%     case 'evaluate'
OutputGenNON = [ UU * coeff];
%    case 'evaluate_derivative'
%        OutputGenNON = [ UU * coeff];
%    case 'evaluate_second_derivative'
%        OutputGenNON = [ UU * coeff];
%end
end
