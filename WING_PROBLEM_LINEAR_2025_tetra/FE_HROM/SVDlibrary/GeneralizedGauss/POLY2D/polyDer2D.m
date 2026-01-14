function [fx,fy] = polyDer2D(p,x,y,n,m)
%POLYDER2D Evaluate derivatives of 2-D polynomial using Horner's method.
%   F = POLYDER2D(P,X,Y,N,M) evaluates the derivatives of 2-D polynomial P at
%   the points specified by X and Y, which must have the same dimensions. The
%   outputs FX and FY will have the same dimensions as X and Y. N and M specify
%   the order of X and Y respectively. Polynomial coefficients are in the
%   following order.
%
%   F(X,Y) = P_1 * X^N * Y^M + P_2 * X^{N-1} * Y^M + ... + P_{N+1} * Y^M + ...
%            P_{N+2} * X^N * Y^{M-1} + P_{N+3} * X^{N-1} * Y^{M-1} + ... + P_{2*(N+1)} * Y^{M-1} + ...
%            ...
%            P_{M*(N+1)+1} * X^N + P_{M*(N+1)+2} * X^{N-1} + ... + P_{(N+1)*(M+1)}
%
% See also: POLYFITN by John D'Errico on MathWorks MATLAB Central FEX
% http://www.mathworks.com/matlabcentral/fileexchange/34765-polyfitn
%
% Mark Mikofski, http://poquitopicante.blogspot.com
% Version 1-0, 2013-04-17

%% LaTex
%
% $$f\left(x,y\right)=p_1 x^n y^m+p_2 x^{\left(n-1\right)} y^m+\ldots+p_{n+1} y^m+\ldots$$
%
% $$p_{n+2}x^ny^{\left(m-1\right)}+p_{n+3}x^{\left(n-1\right)}y^{\left(m-1\right)}+\ldots+p_{2\left(n+1\right)}y^{\left(m-1\right)}+\ldots$$
%
% $$\ldots$$
%
% $$p_{m\left(n+1\right)+1}*x^n+p_{m\left(n+1\right)+2}*x^{\left(n-1\right)}+\ldots+p_{\left(n+1\right)\left(m+1\right)}$$

%% check input args
validateattributes(p,{'numeric'},{'2d','nonempty','real','finite'}, ...
    'polyDer2D','p',1)
validateattributes(x,{'numeric'},{'nonempty','real','finite'}, ...
    'polyDer2D','x',2)
validateattributes(y,{'numeric'},{'nonempty','real','finite'}, ...
    'polyDer2D','y',3)
assert(all(size(x)==size(y)),'polyDer2D:sizeMismatch', ...
    'X and Y must be the same size.')
% use size of p to set n & m
pdims = size(p);
if nargin<4 && all(pdims>1)
    n = pdims(1)-1;
    m = pdims(2)-1;
end
validateattributes(n,{'numeric'},{'scalar','integer','positive','<',10}, ...
    'polyDer2D','n',4)
validateattributes(m,{'numeric'},{'scalar','integer','positive','<',10}, ...
    'polyDer2D','m',5)
if all(pdims>1) 
    assert(pdims(1)==n+1,'polyDer2D:xOrderMismatch', ...
        'The number of x coefficients doesn''t match the order n.')
    assert(pdims(2)==m+1,'polyDer2D:yOrderMismatch', ...
        'The number of y coefficients doesn''t match the order m.')
end
p = p(:);
Np = numel(p);
assert(Np==(n+1)*(m+1),'OrderMismatch', ...
        ['The number of x & y coefficients doesn''t match the orders ', ...
        'n & m.'])
%% fx = df/dx
fx = n*p(1);
for ni = 1:n-1
    fx = fx.*x+(n-ni)*p(1+ni);
end
for mi = 1:m
    mj = (n+1)*mi+1;
    gx = n*p(mj);
    for ni = 1:n-1
        gx = gx.*x+(n-ni)*p(mj+ni);
    end
    fx = fx.*y+gx;
end
%% fy = df/dy
fy = p(1);
for ni = 1:n
    fy = fy.*x+p(1+ni);
end
fy = m*fy;
for mi = 1:m-1
    mj = (n+1)*mi+1;
    gy = p(mj);
    for ni = 1:n
        gy = gy.*x+p(mj+ni);
    end
    fy = fy.*y+(m-mi)*gy;
end