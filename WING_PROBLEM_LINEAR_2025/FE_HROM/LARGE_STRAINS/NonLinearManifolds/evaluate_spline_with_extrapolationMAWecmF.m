function [y,dy,d2y ] = evaluate_spline_with_extrapolationMAWecmF(sp, sp1, sp2, xmin,xmax,INFOextrap,x,VOL)
%--------------------------------------------------------------------------
% This is a modification of evaluate_spline_with_extrapolationFUNv,  (described below) for cases in
% which the function to be approximated is subjected to the following
% constraints:
% 1) sum(sp) = VOL
% 2) sp>= 0
% WE impose this constraints by simply freezing the values at the ends of
% the interval
% Joaquín Alberto Hernández Ortega (JAHO)
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/112_NonLIN_ROM_RBF/10_SAW_ECMlarge.mlx
% 14-Sept-2025, Sunday, 19:43, Balmes 185, Barcelona


% evaluate_spline_with_extrapolationFUNv
%
% PURPOSE
%   Vectorized evaluation of a univariate spline (B-form) and its first and
%   second derivatives at points x, with **quadratic extrapolation** outside
%   [xmin, xmax]. This is a vectorized modification of an existing function
%   (evaluate_spline_with_extrapolation) aimed at reducing loop overhead and
%   minimizing calls to fnval by batching inputs.
%
% WHAT'S NEW (vs. original)
%   - Fully vectorized path: single pass over x using logical masks.
%   - Uses a fast evaluator (spval_cubic_univar) for B-form splines to avoid
%     dispatch overhead; called once per region and derivative order.
%   - Extrapolation computed in bulk with bsxfun/repmat (or implicit
%     expansion in newer MATLAB) to avoid per-point loops.
%
% INPUTS
%   sp      : B-form spline struct for f(x)           (e.g., from spap2)
%   sp1     : B-form spline struct for f'(x)          (typically fnder(sp))
%   sp2     : B-form spline struct for f''(x)         (typically fnder(sp1))
%   xmin    : left boundary of the valid spline interval
%   xmax    : right boundary of the valid spline interval
%   INFOextrap : struct with Taylor data at boundaries (one column per output)
%                .XMIN.s0  = f(xmin)   size: [nfun × 1]
%                .XMIN.s1  = f'(xmin)  size: [nfun × 1]
%                .XMIN.s2  = f''(xmin) size: [nfun × 1]
%                .XMAX.s0  = f(xmax)   size: [nfun × 1]
%                .XMAX.s1  = f'(xmax)  size: [nfun × 1]
%                .XMAX.s2  = f''(xmax) size: [nfun × 1]
%   x       : row/column vector or array of evaluation points
%
% OUTPUTS
%   y   : f(x)   evaluated at all points in x (same shape as x, per row-block)
%   dy  : f'(x)  evaluated at all points in x
%   d2y : f''(x) evaluated at all points in x
%
% ASSUMPTIONS / CONTRACT
%   - Univariate, scalar- or multi-output spline in B-form. Number of output
%     functions nfun = size(sp.coefs,1). The same nfun applies to sp1/sp2.
%   - sp1 = fnder(sp), sp2 = fnder(sp1) (consistent derivatives).
%   - xmin < xmax, and x may lie anywhere on the real line.
%   - Extrapolation is **quadratic**:  f(xb+dx) ≈ s0 + s1*dx + 0.5*s2*dx^2,
%     with (s0,s1,s2) taken from INFOextrap at xb ∈ {xmin, xmax}.
%
% NUMERICAL DETAILS
%   - In-domain points use spval_cubic_univar (a lightweight B-form evaluator).
%   - Left/right out-of-domain points use bulk Taylor evaluation with bsxfun
%     (compatible with older MATLAB). In R2016b+, bsxfun/repmat can be
%     replaced by implicit expansion: s1.*dx, s2.*dx, s2.*(dx.^2), etc.
%
% PERFORMANCE NOTES
%   - Single allocation of y,dy,d2y with vectorized fills per region.
%   - Three evaluator calls for in-domain region (f,f',f'') amortize overhead.
%   - Extrapolation is memory-bound; chunk x if extremely large to cap RAM.
%
% EXAMPLE
%   sp  = spap2(augknt(knots,4), 4, xdata, ydata);
%   sp1 = fnder(sp); sp2 = fnder(sp1);
%   INFOextrap.XMIN.s0 = fnval(sp, xmin);   INFOextrap.XMAX.s0 = fnval(sp, xmax);
%   INFOextrap.XMIN.s1 = fnval(sp1,xmin);   INFOextrap.XMAX.s1 = fnval(sp1,xmax);
%   INFOextrap.XMIN.s2 = fnval(sp2,xmin);   INFOextrap.XMAX.s2 = fnval(sp2,xmax);
%   [y,dy,d2y] = evaluate_spline_with_extrapolationFUNv(sp,sp1,sp2,xmin,xmax,INFOextrap,x);
%
% LIMITATIONS
%   - No left-continuity option (right-continuous spline convention).
%   - Assumes derivative splines (sp1, sp2) are provided and consistent.
%   - No input validation beyond basic shape handling (for speed).
%
% AUTHOR / HISTORY
%   Vectorized version by Joaquín A. Hernández, 11-aG-2025 (cARTAGENA).
%   Based on a prior non-vectorized routine "evaluate_spline_with_extrapolation".
%--------------------------------------------------------------------------

% --- Shape normalization and preallocation (vectorized output buffers)

%--------------------------------------------------------------------------
if nargin == 0
    load('tmp.mat')
    x = [-20,20] ;
end
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
%
x = x(:)';  % Ensure row vector
nfun = size(sp.coefs,1) ;
y   = zeros(nfun,length(x));
dy  = y ;
d2y = y ;

% % Define domain
% interval_t= fnbrk(sp, 'interval');
% xmin = interval_t(1);
% xmax = interval_t(2);

% Logical masks for domain regions
in  = (x >= xmin) & (x <= xmax);
left  = x < xmin;
right = x > xmax;

INDICES = 1:length(x) ;
in = INDICES(in) ;
left = INDICES(left) ;
right = INDICES(right) ;


% Evaluate in-domain using fnval
% y(in)   = spval(sp,   x(in));
%  dy(in)  = spval(sp1,  x(in));
%  d2y(in) = spval(sp2,  x(in));

if ~isempty(in)
    y(:,in)   = spval_cubic_univar(sp,   x(in));
    dy(:,in)  = spval_cubic_univar(sp1,  x(in));
    d2y(:,in) = spval_cubic_univar(sp2,  x(in));
end


% Left extrapolation (quadratic)
if ~isempty(left)
    
    dx = x(left) - xmin;
    s0 = repmat(INFOextrap.XMIN.s0,1,length(dx));
%     s1 = repmat(INFOextrap.XMIN.s1,1,length(dx));
%     s2 = repmat(INFOextrap.XMIN.s2,1,length(dx));
    
%     s1_dx = bsxfun(@times,s1',dx')' ;
%     s2_dx = bsxfun(@times,s2',dx')' ;
%     s2_dx2 = bsxfun(@times,s2',(dx.^2'))' ;
    
    % y(:,left)   = s0 + s1 * dx + 0.5 * s2 * dx.^2;
    % dy(:,left)  = s1 + s2 * dx;
    % d2y(:,left) = s2;
    
    y(:,left)   = s0  ; %+ s1_dx + 0.5 * s2_dx2;
    dy(:,left)  = zeros(size(s0)) ;% + s2_dx;
    d2y(:,left) = zeros(size(s0));; %
    
    % Project each extrapolated column onto the simplex (sum=VOL, >=0)
  %  [y,dy,d2y] = ProjectActiveSet_ECM(y,left,VOL,dy,d2y) ; 
    
end

% Right extrapolation (quadratic)
if any(right)
    dx = x(right) - xmax;
    s0 = repmat(INFOextrap.XMAX.s0,1,length(dx));
    s1 = repmat(INFOextrap.XMAX.s1,1,length(dx));
    s2 = repmat(INFOextrap.XMAX.s2,1,length(dx));
    
 %   s1_dx = bsxfun(@times,s1',dx')' ;
 %   s2_dx = bsxfun(@times,s2',dx')' ;
 %   s2_dx2 = bsxfun(@times,s2',(dx.^2'))' ;
    
    
    y(:,right)   = s0 ; %+ s1_dx + 0.5 * s2_dx2;
    dy(:,right)  =  zeros(size(s0)) ; %s1 + s2_dx;
    d2y(:,right) = zeros(size(s0)) ; % s2;
    
    
 %   [y,dy,d2y] = ProjectActiveSet_ECM(y,right,VOL,dy,d2y) ; 
end

