function c = multinomial(n,varargin)
% MULTINOMIAL Multinomial coefficients
%
%   MULTINOMIAL(N, K1, K2, ..., Km) where N and Ki are floating point
%   arrays of non-negative integers satisfying N = K1 + K2 + ... + Km, 
%   returns the multinomial  coefficient   N!/( K1!* K2! ... *Km!).
%
%   MULTINOMIAL(N, [K1 K2 ... Km]) when Ki's are all scalar, is the 
%   same as MULTINOMIAL(N, K1, K2, ..., Km) and runs faster.
%
%   Non-integer input arguments are pre-rounded by FLOOR function.
%
% EXAMPLES:
%    multinomial(8, 2, 6) returns  28 
%    binomial(8, 2) returns  28
% 
%    multinomial(8, 2, 3, 3)  returns  560
%    multinomial(8, [2, 3, 3])  returns  560
%
%    multinomial([8 10], 2, [6 8]) returns  [28  45]

% Mukhtar Ullah
% November 1, 2004
% mukhtar.ullah@informatik.uni-rostock.de

nIn = nargin;
error(nargchk(2, nIn, nIn))

if ~isreal(n) || ~isfloat(n) || any(n(:)<0)
    error('Inputs must be floating point arrays of non-negative reals')
end

arg2 = varargin; 
dim = 2;

if nIn < 3                         
    k = arg2{1}(:).'; 
    if isscalar(k)
        error('In case of two arguments, the 2nd cannot be scalar')
    end    
else
    [arg2{:},sizk] = sclrexpnd(arg2{:});
    if sizk == 1
        k = [arg2{:}];        
    else
        if ~isscalar(n) && ~isequal(sizk,size(n))
            error('Non-scalar arguments must have the same size')
        end
        dim = numel(sizk) + 1; 
        k = cat(dim,arg2{:});              
    end    
end

if ~isreal(k) || ~isfloat(k) || any(k(:)<0)
    error('Inputs must be floating point arrays of non-negative reals')
end

n = floor(n);
k = floor(k);

if any(sum(k,dim)~=n)
    error('Inputs must satisfy N = K1 + K2 ... + Km ')
end

c = floor(exp(gammaln(n+1) - sum(gammaln(k+1),dim)) + .5); 