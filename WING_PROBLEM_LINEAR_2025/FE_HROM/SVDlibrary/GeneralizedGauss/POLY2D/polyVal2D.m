function f = polyVal2D(p,x,y,n,m)
%POLYVAL2D Evaluate a 2-D polynomial using Horner's method.
%   F = POLYVAL2D(P,X,Y) evaluates the 2-D polynomial P at the points 
%   specified by X and Y, which must have the same dimensions. The output F
%   will have the same dimensions as X and  Y. N and M specify the order of X
%   and Y respectively. Polynomial coefficients are in the following order.
%
%   F(X,Y) = P_1 * X^N * Y^M + P_2 * X^{N-1} * Y^M + ... + P_{N+1} * Y^M + ...
%            P_{N+2} * X^N * Y^{M-1} + P_{N+3} * X^{N-1} * Y^{M-1} + ... + P_{2*(N+1)} * Y^{M-1} + ...
%            ...
%            P_{M*(N+1)+1} * X^N + P_{M*(N+1)+2} * X^{N-1} + ... + P_{(N+1)*(M+1)}
%
%
% See also: POLYFITN by John D'Errico on MathWorks MATLAB Central FEX
% http://www.mathworks.com/matlabcentral/fileexchange/34765-polyfitn
%
% Mark Mikofski, http://poquitopicante.blogspot.com
% Version 1-0, 2013-03-13

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
    'polyVal2D','p',1)
validateattributes(x,{'numeric'},{'nonempty','real','finite'}, ...
    'polyVal2D','x',2)
validateattributes(y,{'numeric'},{'nonempty','real','finite'}, ...
    'polyVal2D','y',3)
assert(all(size(x)==size(y)),'polyVal2D:sizeMismatch', ...
    'X and Y must be the same size.')
% use size of p to set n & m
pdims = size(p);
if nargin<4 && all(pdims>1)
    n = pdims(1)-1;
    m = pdims(2)-1;
end
validateattributes(n,{'numeric'},{'scalar','integer','positive','<',10}, ...
    'polyVal2D','n',4)
validateattributes(m,{'numeric'},{'scalar','integer','positive','<',10}, ...
    'polyVal2D','m',5)
if all(pdims>1) 
    assert(pdims(1)==n+1,'polyVal2D:xOrderMismatch', ...
        'The number of x coefficients doesn''t match the order n.')
    assert(pdims(2)==m+1,'polyVal2D:yOrderMismatch', ...
        'The number of y coefficients doesn''t match the order m.')
end
%% evaluate polynomial P
f = p(1);
for ni = 1:n
    f = f.*x+p(1+ni);
end
for mi = 1:m
    mj = (n+1)*mi+1;
    g = p(mj);
    for ni = 1:n
        g = g.*x+p(mj+ni);
    end
    f = f.*y+g;
end
