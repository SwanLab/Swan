clc,clear,close all
E = 210; nu = 0.3;
Gc0 = 5e-3; sigma = 1.5;
%% l0 constant
% sigmach = 1.5;
% Gc = Gc0*(sigma^2/sigmach^2);
% l0 = 0.1;

%% l0 variable
Gc = Gc0;
lch = (2*E*Gc)/(sigma^2);

C11 = E/((1+nu)*(1-nu));
k  = E./(2.*(1-nu));
mu = E./(2.*(1+nu));
slope = (- k - k^2/mu - mu - (2*mu^2 + k)/k)/C11;
lhs = -2*(3/8)*(Gc/slope)*(E/sigma^2);

l0 = lhs;

%% Inputs
mat.E = E;
mat.sigma = sigma;
mat.Gc = Gc;
mat.l0 = l0;

nCoeffs = [9 5 9 5 9 5];
data = 'SquareArea';
output = 'SquareAreaLowDegreeHashimStrikman';
fittingLessDenominator(data,mat,nCoeffs,output)