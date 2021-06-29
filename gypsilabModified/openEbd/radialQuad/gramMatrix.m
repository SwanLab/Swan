function [ A,T,a,timer ] = gramMatrix(a,rho)
% A  = gramMatrix(a,rho)
% Inputs : 
% - a : (scalar / array of size 2). If a = [a1,a2], the matrix contains the
% scalar products on the ring {a1 < r < a2}. If a is scalar, a2 is taken as
% 1. 
% - rho : array of positive reals. 
% Output : 
% Returns the square matrix defined A(p,q) = (e_p | e_q), more precisely,
% A_{p,q} = 2\pi \int_{a(1)}^{a(2)} r (e_p'(r) e_q'(r))dr
% where e_p(r) = C_p J_0(rho(p) r) are eigenfunctions of the Laplace
% operator, and C_p is the normalization constant (see function Cp)
% 'rho' is a positive vector of frequencies. In the classical setting, it
% contains roots of the function J_0 but the formulas are also valid for
% arbitrary rho.
% the parameter 'a' can be a scalar (then a(2) is replaced by 1).
% [ A,T ] = gramMatrix(a,rho) returns the cholesky factorization of A, that
% is T is such that A = T'*T.
% [ A,T,timer] = gramMatrix(a,rho) returns the time it took to compute the
% matrix A and its cholesky factorization.

tic
rho = rho(:);
P = length(rho);
if length(a) == 2
    b = a(2);
    a = a(1);
else
    b = 1;
end
if b == 1
    if a > 8/P
        warning(...
            ['Parameer ''a'' is likely too high, '...
            'A will be ill-conditioned, consider decreasing']...
            )
        % Valid only if rho is ordered as usual.
    end
end

% Extra diagonal values
C = Cp(rho);
ek = @(r)(C.*besselj(0,rho*r));
el = @(r)(ek(r)');
dek = @(r)(C.*rho.*besselj(1,rho*r));
del = @(r)(dek(r)');
rhok = rho;
rhol = rho';
U = (rhok*(1./rhol)).^2;
V = 1./U;

F = @(r)( r* ((dek(r)*el(r))./(U-1) + (ek(r)*del(r))./(V-1)) );
A = 2*pi*(F(b) - F(a));

% Diagonal values

D = @(r)(1/2*C.^2.*rho.^2.*r^2.*(besselj(1,rho*r).^2 - besselj(0,rho*r).*besselj(2,rho*r)));
A(1:P+1:end) = 2*pi*(D(b) - D(a));


if nargout >=2
    try
        T = chol(A);
    catch
        warning('Condition number too high, restarting with a = %d',a)
        a = 3/4*a;
        timer1 = toc;
        [A,T,a,timer2] = gramMatrix(a,rho);
        timer = timer1 + timer2;
        return
    end
end

timer = toc;

end

