function [OUTPUT ] = evaluate_spline_with_extrapolationECMw(sp, sp1, sp2, x)
% This is an adaptation of evaluate_spline_with_extrapolation (described below) for cases in
% which the function to be approximated is subjected to the following
% constraints:
% 1) sum(sp) = VOL
% 2) sp>= 0
% Joaquín Alberto Hernández Ortega (JAHO)
% 14-Sept-2025, Sunday, 16:00, Balmes 185, Barcelona

%--------------------------------------------------------------------------
% function [OUTPUT] = evaluate_spline_with_extrapolation(sp, sp1, sp2, x).
%
% PURPOSE:
%   Evaluates the value, first derivative, and second derivative of a spline
%   function at specified points `x`. When `x` lies outside the spline's domain,
%   the function performs **quadratic extrapolation** using Taylor expansion
%   at the nearest boundary (left or right).
%
% INPUTS:
%   - sp   : spline function in B-form, e.g., obtained from `spap2`
%   - sp1  : first derivative of `sp`, typically `fnder(sp)`
%   - sp2  : second derivative of `sp`, typically `fnder(sp1)`
%   - x    : vector of evaluation points (can be inside or outside the domain)
%
% OUTPUT:
%   - OUTPUT : structure containing:
%       • OUTPUT.VALUE : values of the spline at x
%       • OUTPUT.DER1  : first derivatives at x
%       • OUTPUT.DER2  : second derivatives at x
%
% NOTES:
%   - Inside the spline domain, the function uses `fnval` for evaluation.
%   - Outside the domain, it extrapolates using a second-order Taylor
%     expansion around the closest boundary point (xmin or xmax).
%
% EXAMPLE USAGE:
%   sp = spap2(knots, degree, xdata, ydata);
%   sp1 = fnder(sp);
%   sp2 = fnder(sp1);
%   x_eval = linspace(...);
%   OUT = evaluate_spline_with_extrapolation(sp, sp1, sp2, x_eval);
%
% AUTHOR:
%   Joaquín A. Hernández, 12-July-2025, Origins, Calvet Barcelona
%--------------------------------------------------------------------------

%
% Inputs:
%   sp   - spline in B-form (from spap2)
%   sp1  - first derivative (fnder(sp))
%   sp2  - second derivative (fnder(sp1))
%   x    - evaluation points (vector)
%
% Outputs:
%   y    - spline value at x (with extrapolation)
%   dy   - first derivative at x
%   d2y  - second derivative at x

if nargin == 0
    load('tmp.mat')
end

x = x(:)';  % Ensure row vector
y   = zeros(size(x));
dy  = zeros(size(x));
d2y = zeros(size(x));

% Define domain
interval_t= fnbrk(sp, 'interval');
xmin = interval_t(1);
xmax = interval_t(2);

% Logical masks for domain regions
in  = (x >= xmin) & (x <= xmax);
left  = x < xmin;
right = x > xmax;

% Evaluate in-domain using fnval
y(in)   = fnval(sp,   x(in));
dy(in)  = fnval(sp1,  x(in));
d2y(in) = fnval(sp2,  x(in));

% Left extrapolation (quadratic)
if any(left)
    dx = x(left) - xmin;
    s0 = fnval(sp, xmin);
    s1 = fnval(sp1, xmin);
    s2 = fnval(sp2, xmin);
    
    y(left)   = s0 + s1 * dx + 0.5 * s2 * dx.^2;
    dy(left)  = s1 + s2 * dx;
    d2y(left) = s2;
end

% Right extrapolation (quadratic)
if any(right)
    dx = x(right) - xmax;
    s0 = fnval(sp, xmax);
    s1 = fnval(sp1, xmax);
    s2 = fnval(sp2, xmax);
    
    y(right)   = s0 + s1 * dx + 0.5 * s2 * dx.^2;
    dy(right)  = s1 + s2 * dx;
    d2y(right) = s2;
end
OUTPUT.VALUE = y ;
OUTPUT.DER1 = dy ;
OUTPUT.DER2 = d2y ;

end
