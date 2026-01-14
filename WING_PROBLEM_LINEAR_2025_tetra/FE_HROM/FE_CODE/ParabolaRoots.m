function roots = ParabolaRoots(varargin)
% function roots = ParabolaRoots(P)
%
% Find roots of second order polynomials
%
% INPUT
%   P: (n x 2) array, each row corresponds to coefficients of each
%   polynomial, P(:,1)*x^2 + P(:,2)*x + P(:,3)
% OUTPUT
%   roots: (n x 2) array, each row correspond to the roots of P
%
% To adjust the parameter below which the the discriminant is considerered
% as nil, use
%   ParabolaRoots(P, tol)
% Adjusting tol is useful to avoid the real roots become complex due to
% numerical accuracy. The default TOL is 0
%
% See also: roots, CardanRoots, eig2
%
% Author: Bruno Luong <brunoluong@yahoo.com>
% History:
%     Original 27-May-2010

% Adjustable parameter
tol = 0;

if nargin<3
    P = varargin{1};
    a = P(:,1);
    b = P(:,2);
    c = P(:,3);
    if nargin>=2
        tol = varargin{2};
    end
else
    [a b c] = deal(varargin{1:3});
    if nargin>=4
        tol = varargin{2};
    end
end

if ~isequal(a,1)
    b = b./a;
    c = c./a;
end

b = 0.5*b;

delta = b.^2 - c;
delta(abs(delta)<tol) = 0;
sqrtdelta = sqrt(delta);

roots = [sqrtdelta -sqrtdelta];
roots = bsxfun(@minus, roots, b);

end

