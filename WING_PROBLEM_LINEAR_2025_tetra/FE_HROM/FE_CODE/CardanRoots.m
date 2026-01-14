function roots = CardanRoots(varargin)
% function roots = CardanRoots(P)
%
% Find roots of third order polynomials using Cardan's formula
%
% INPUT
%   P: (n x 4) array, each row corresponds to coefficients of each
%   polynomial, P(:,1)*x^3 + P(:,2)*x^2 + P(:,3)*x + P(:,4)
% OUTPUT
%   roots: (n x 3) array, each row correspond to the roots of P
%
% To adjust the parameter below which the the discriminant is considerered
% as nil, use
%   CardanRoots(P, tol)
% Adjusting tol is useful to avoid the real roots become complex due to
% numerical accuracy. The default TOL is 0
%
% http://www.sosmath.com/algebra/factor/fac11/fac11.html
% http://mathforum.org/dr.math/faq/faq.cubic.equations.html
%
% See also: roots, ParabolaRoots, eig3
%
% Author: Bruno Luong <brunoluong@yahoo.com>
% History:
%     Original 20-May-2010

% Adjustable parameter
tol = 0;

if nargin<4
    P = varargin{1};
    a = P(:,1);
    b = P(:,2);
    c = P(:,3);
    d = P(:,4);
    if nargin>=2
        tol = varargin{2};
    end
else
    [a b c d] = deal(varargin{1:4});
    if nargin>=5
        tol = varargin{5};
    end
end

if ~isequal(a,1)
    b = b./a;
    c = c./a;
    d = d./a;
end

b2 = b.^2;

p = -b2/3 + c;
q = ((2/27)*b2-(1/3)*c).*b + d;
delta = q.^2 + (4/27)*p.^3;

% Three cases of discriminant sign
iscmplx = imag(p) | imag(q);
notcmplx = ~iscmplx;
deltanull = notcmplx & abs(delta)<tol; % = 0
deltaneg = notcmplx & delta<0;
deltapos = notcmplx & ~(deltanull | deltaneg);

n = size(delta,1);
roots = zeros(n, 3, class(delta));

%%
if any(deltanull)
    idx = find(deltanull);
    roots(idx,:) = CardanNull(p(idx), q(idx));
end

%%
if any(deltaneg)
    idx = find(deltaneg);
    roots(idx,:) = CardanNeg(q(idx), delta(idx));
end

%%
if any(deltapos)
    idx = find(deltapos);
    roots(idx,:) = CardanPos(q(idx), delta(idx));
end

%%
if any(iscmplx)
    idx = find(iscmplx);
    roots(idx,:) = CardanCmplx(p(idx), q(idx), delta(idx));
end

%% Common for all cases
roots = bsxfun(@minus, roots, b/3);

end

%%
function roots = CardanNull(p, q)

S1 = 3*q./p;
S2 = -0.5*S1;
roots = [S1 S2 S2];  % double real solutions

end % CardanNull

%%
function roots = CardanNeg(q, delta)

alfa = -q;
beta = sqrt(-delta);
r2 = alfa.^2-delta;
rho = (4^(1/3))*exp(log(r2)/6);
theta = atan2(beta,alfa)/3;
alfa = rho.*cos(theta);
beta = rho.*sin(theta);
S1 = alfa;
x = (-0.5)*alfa;
y = (sqrt(3)/2)*beta;
S2 = x-y;
S3 = x+y;
roots = [S1 S2 S3];

end % CardanNeg

%%
function roots = CardanPos(q, delta)

sqrtdelta = sqrt(delta);
u3 = (-q+sqrtdelta)/2;
v3 = (-q-sqrtdelta)/2;

% Cubic roots of u3 and v3
u = sign(u3).*exp(log(abs(u3))/3);
v = sign(v3).*exp(log(abs(v3))/3);

S1 = u+v;
% Complex solutions
j = complex(-0.5,sqrt(3)/2);
j2 = complex(-0.5,-sqrt(3)/2);
S2 = j*u+j2*v;
S3 = conj(S2);
roots = [S1 S2 S3];

end % CardanPos

%%
function roots = CardanCmplx(p, q, delta)

sqrtdelta = sqrt(delta);
u3 = (-q+sqrtdelta)/2;
v3 = (-q-sqrtdelta)/2;

% we need u*v = -p/3
p = (-1/3)*p;
iu = abs(u3)>abs(v3);
u = zeros(size(u3),class(u3));
v = zeros(size(v3),class(v3));
u(iu) = exp(log(u3(iu))/3);
v(iu) = p(iu)./u(iu);
v(~iu) = exp(log(v3(~iu))/3);
u(~iu) = p(~iu)./v(~iu);

S1 = u+v;

j = complex(-0.5,sqrt(3)/2);
j2 = complex(-0.5,-sqrt(3)/2);
S2 = j*u+j2*v;
S3 = j2*u+j*v;
roots = [S1 S2 S3];

end % CardanCmplx

