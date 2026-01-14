function p = polyFit2D(f,x,y,n,m)
%POLYFIT2D Fit data with a 2-D polynomial.
%   P = POLYFIT2D(F,X,Y,N,M) determines the 2-D polynomial P that best fit the
%   data F(X,Y) in the least squares sense. Polynomial coefficients are in the
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
validateattributes(f,{'numeric'},{'nonempty','real','finite'}, ...
    'polyFit2D','f',1)
validateattributes(x,{'numeric'},{'nonempty','real','finite'}, ...
    'polyFit2D','x',2)
validateattributes(y,{'numeric'},{'nonempty','real','finite'}, ...
    'polyFit2D','y',3)
assert(all(size(x)==size(y)),'polyFit2D:sizeMismatch', ...
    'X and Y must be the same size.')
assert(all(size(x)==size(f)),'polyFit2D:sizeMismatch', ...
    'X and F must be the same size.')
validateattributes(n,{'numeric'},{'scalar','integer','positive','<',10}, ...
    'polyFit2D','n',4)
validateattributes(m,{'numeric'},{'scalar','integer','positive','<',10}, ...
    'polyFit2D','m',5)
Npoints = numel(f);
assert(all((n+1)*(m+1)<Npoints),'polyFit2D:degreeTooBig', ...
    'Degree must be smaller than number of data poiints.')
%% fit data
f = f(:);
x = x(:)*ones(1,n+1);
y = y(:)*ones(1,(n+1)*(m+1));
z = ones(Npoints,(n+1)*(m+1));
for ni = 1:n
    z(:,1:ni) = z(:,1:ni).*x(:,1:ni);
end
for mi = 1:m
    mj = (n+1)*mi;
    z(:,1:mj) = z(:,1:mj).*y(:,1:mj);
    for ni = 1:n
        z(:,mj+1:mj+ni) = z(:,mj+1:mj+ni).*x(:,1:ni);
    end
end
p = z\f;