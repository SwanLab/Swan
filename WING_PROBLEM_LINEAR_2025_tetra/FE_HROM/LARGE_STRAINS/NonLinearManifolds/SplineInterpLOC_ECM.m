function [NonLinear_q, NonLinear_q_der,nREDcoor] = SplineInterpLOC_ECM(DATA_interp, qINF, qSUP, UU, SS)
% SplineInterpLOC - Constructs interpolated spline models from qSUP vs qINF,
% and returns two nonlinear evaluators:
%   NonLinear_q(x): evaluates [x; UU * (SS .* spline(x))]
%   NonLinear_q_der(x): evaluates [1; UU * (SS .* spline'(x))]. 
%
% Inputs:
%   DATA_interp - struct with field NSAMPLES: number of points for interpolation
%   qINF        - 1xN vector of input (abscissas)
%   qSUP        - mxN matrix of output values (each row = function)
%   UU          - matrix (typically POD modes or basis) used to project spline values
%   SS          - vector of scaling factors for each output mode
%
% Outputs:
%   NonLinear_q       - function handle: evaluates spline projections at x
%   NonLinear_q_der   - function handle: evaluates spline derivative projections at x

if nargin == 0
    load('tmp1.mat')  % Fallback for debugging/demo
end

nREDcoor = 1;

% Validate inputs
assert(isvector(qINF), 'qINF must be a vector');
[m, N] = size(qSUP);
assert(length(qINF) == N, 'qINF and qSUP must have compatible dimensions');

% Uniform sampling of n points from qINF and qSUP
% n = DATA_interp.NSAMPLES;
% qmin = qINF(1) ; qmax = qINF(end);
% qINF_red_approx = linspace(qmin, qmax, n) ;
% [ idx] = knnsearch(qINF(:),qINF_red_approx(:)) ;
% idx = unique(idx) ;
qINF_red = qINF ;
qSUP_red = qSUP  ;
%
%
%     idx = round(linspace(1, N, n));
%     qINF_red = qINF(idx);
%     qSUP_red = qSUP(:, idx);

% Initialize spline model storage
SplineModel = struct('spline', {}, 'derivative', {}, ...
    'evaluate', {}, 'evaluate_derivative', {});

% Plot setup
figure(150);
hold on;
colors = lines(m);
xlabel('Amst');
ylabel('Aslv');
title('Spline Fitting of Aslv vs Amst');
hold off;
% 
% figure(20);
% hold on
% 
% xlabel('qINF');
% ylabel('d qSUP/dq');
% %      hold off;
%
%       sp = csaps(qINF_red, y, p);
DATA_interp = DefaultField(DATA_interp,'RegularizationTermSplines',[]);

% Build splines for each output function
for i = 1:m
    y = qSUP_red(i, :);
    
    % Create spline and its derivative
    
    
    if isempty(DATA_interp.RegularizationTermSplines)
        sp = spline(qINF_red, y);
    else
        sp = csaps(qINF_red, y, DATA_interp.RegularizationTermSplines);
    end
    sp_der = fnder(sp);
    
    % Store evaluators in structure
    SplineModel(i).spline = sp;
    SplineModel(i).derivative = sp_der;
    SplineModel(i).evaluate = @(x) ppval(sp, x);
    SplineModel(i).evaluate_derivative = @(x) ppval(sp_der, x);
    
    % Plot original and fitted curve
    figure(150)
    hold on
    yy = SplineModel(i).evaluate(qINF);
    plot(qINF, qSUP(i,:),'.','Color', colors(i,:), ...
        'DisplayName', ['\lambda_{', num2str(i),'} = ', num2str(SS(i))]);
    
    plot(qINF, yy, '--', 'Color', colors(i,:), ...
        'DisplayName', ['Spline ', num2str(i)]);
    
%     %
%     plot(qINF_red,qSUP_red(i,:),'Marker','x', 'Color', colors(i,:), ...
%         'DisplayName', ['Spline P. ', num2str(SS(i))],'LineStyle','none');
%     
    
%     figure(16)
%     hold on
%     dydx = SplineModel(i).evaluate_derivative(qINF);
%     plot(qINF, dydx, '-', ...
%         'DisplayName', ['d\lambda_{', num2str(i), '}/dq']);
    
    
end
% figure(160)
% legend show;
figure(150)
legend show;


% Generate function handles for evaluating spline-based reduced outputs
[NonLinear_q, NonLinear_q_der] = make_spline_evaluator(SplineModel, UU, SS);
end

function [NonLinear_q, NonLinear_q_der] = make_spline_evaluator(SplineModel, UU, SS)
% make_spline_evaluator - Generates two function handles:
%   NonLinear_q(x)       - evaluates [x; UU * (SS .* spline(x))]
%   NonLinear_q_der(x)   - evaluates [1; UU * (SS .* spline'(x))]
%
% Inputs:
%   SplineModel - array of spline models with .evaluate and .evaluate_derivative
%   UU          - projection matrix (e.g., from POD or Galerkin)
%   SS          - scaling factors for each function
%
% Outputs:
%   NonLinear_q       - nonlinear evaluation function (value)
%   NonLinear_q_der   - nonlinear evaluation function (derivative)

NonLinear_q = @(x) fast_eval2(SplineModel, x, 'evaluate', UU, SS);
NonLinear_q_der = @(x) fast_eval2(SplineModel, x, 'evaluate_derivative', UU, SS);
end


function OutputGenNON = fast_eval2(SplineModel, x, fieldnameLABEL, UU, SS)
% fast_eval - Evaluates a set of spline functions or their derivatives at x,
% applies scaling SS(i), and projects the result via UU
%
% Inputs:
%   SplineModel     - structure array with .evaluate / .evaluate_derivative
%   x               - evaluation point(s), scalar or vector
%   fieldnameLABEL  - 'evaluate' or 'evaluate_derivative'
%   UU              - projection matrix (e.g., POD basis)
%   SS              - scaling factors (1 per spline function)
%
% Output:
%   OutputGenNON    - vector:
%       if fieldnameLABEL = 'evaluate'           → [x; UU * scaled_spline_values]
%       if fieldnameLABEL = 'evaluate_derivative'→ [1; UU * scaled_derivatives]

x = x(:).';  % Ensure x is a row vector
m = numel(SplineModel);
n = numel(x);

coeff = zeros(m, n);  % Store scaled spline or derivative values

for i = 1:m
    % Evaluate either spline or its derivative, and scale by SS(i)
    coeff(i, :) = SS(i) * SplineModel(i).(fieldnameLABEL)(x);
end

switch fieldnameLABEL
    case 'evaluate'
        OutputGenNON = [x; UU * coeff];  % Return [x; projected value]
    case 'evaluate_derivative'
        OutputGenNON = [1; UU * coeff];  % Return [1; projected derivative]
end
end
