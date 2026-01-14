function [fnew,dfnew,a,da,xSAMPL] = InterpolatePolyOver(x,F,xnew,m,N)

% 1D Polynomial  fitting of order m using the N closest points to xnew (within
% x). F are the values of the function
% JAHO, March-29th,2020

% STEP 1
% ------
% Detect the set of N points closest to xnew
if length(x) < N ; error('The number of sampling points N is too high') ; end

% STEP 2
% ------
% Identify the N points that  are  closest to xnew
dx = abs(x- xnew);
[~,IND] = sort(dx);
IND = IND(1:N) ;
xSAMPL = x(IND)' ;
[xSAMPL,III] = sort(xSAMPL) ;
fSAMPL = F(IND)' ;
fSAMPL = fSAMPL(III);

% STEP 3
%-----
% Determine the coefficients of the polynomial
C = ones(N,m+1) ;
for i = 2:m+1
    C(:,i) = xSAMPL.^(i-1) ;
end
a = C\fSAMPL ;

a =a((m+1):-1:1) ;  % Coefficients in descending order

% Derivative
nder = length(a)-1 ;
coeff = m:-1:1;
da = a(1:m).*coeff';

fnew = polyval(a,xnew) ;
dfnew = polyval(da,xnew) ;