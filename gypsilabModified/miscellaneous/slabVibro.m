function [R,T] = slabVibro(f, rho, c, e)
% Copyright (c) 2018-2019, Marc Bakry, Ecole Polytechnique       
% GNU General Public License v3.0. 
% Computation of the coeff. of reflection / transmission for an acaustic wave
% crossing three environments with the parameters rho_ {1, 2, 3} and 
% c_ {1,2, 3} for a set of frequencies f
% e.g. ~/nonRegressionTest/vibroAcoustic/transmission_1d_calcul
% To be commented...
assert(length(rho) == 3 && size(c, 1) == length(f) && size(c, 2) == 3)
omega = 2*pi*f;
R = zeros(size(f)); T = zeros(size(f));
muk = (ones(length(f), 1)*rho) .* c;
for i=1:length(f)
    tmp = muk(i, :)*omega(i);
    k   = omega(i)./c(i, :);
    b   = [-1; tmp(1); 0; 0];
    A   = [1 -1 -1 0; ...
        tmp(1) -tmp(2) tmp(2) 0; ...
        0 exp(-1i*k(2)*e) exp(1i*k(2)*e) -exp(1i*k(3)*e); ...
        0 tmp(2)*exp(-1i*k(2)*e) ...
        -tmp(2)*exp(1i*k(2)*e) tmp(3)*exp(1i*k(3)*e)];
    x = A\b;
    R(i) = x(1);
    T(i) = x(4);
end
end