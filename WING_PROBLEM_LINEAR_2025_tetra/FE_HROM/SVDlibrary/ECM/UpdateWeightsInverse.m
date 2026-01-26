function [Bast x] = UpdateWeightsInverse(A,Aast,a,xold,r)
%--------------------------------------------------------------------------
% function [Bast, x] = UpdateWeightsInverse(A, Aast, a, xold, r)
%
% PURPOSE:
%   Efficiently updates the inverse of a Gram matrix and the least-squares 
%   solution when a new column `a` is appended to matrix `A`. This is done 
%   without recomputing the full inverse, using a low-rank update formula.
%
%   The update is particularly useful in iterative methods such as DEIM/ECM 
%   or greedy basis selection, where the matrix grows by one column at each step.
%
% INPUTS:
%   - A     : Original matrix [P x k], where columns are basis vectors
%   - Aast  : Precomputed inverse of (AᵗA), i.e., Aast = inv(A'*A)
%   - a     : New column vector [P x 1] to be added to A
%   - xold  : Previous least-squares solution, xold = Aast * (Aᵗ * b)
%   - r     : Residual vector, r = b - A * xold
%
% OUTPUTS:
%   - Bast  : Updated inverse of ( [A a]ᵗ * [A a] ), computed via Sherman–Morrison–Woodbury
%   - x     : Updated least-squares solution for [A a], i.e., x = inv([A a]' * [A a]) * [A a]' * b
%
% INTERNAL COMPUTATIONS:
%   - Uses block inverse formula for symmetric positive-definite matrices:
%       Bast = [ Aast + (d*dᵗ)/s , -d/s ;
%                -dᵗ/s           , 1/s ]
%     where:
%       c = Aᵗa, d = Aast * c, s = aᵗa - cᵗ * d
%   - Updated solution x is then built from xold and v = (aᵗ * r) / s
%
% EXAMPLE (from test case inside nargin == 0 block):
%   P = 100; k = 49;
%   A = randn(P,k); a = randn(P,1);
%   Aast = inv(A'*A);
%   b = randn(P,1);
%   xold = Aast*(A'*b); r = b - A*xold;
%   [Bast, x] = UpdateWeightsInverse(A, Aast, a, xold, r);
%
% AUTHOR:
%   Joaquín A. Hernández, UPC-CIMNE, May 2025.
%--------------------------------------------------------------------------

if nargin == 0
%     format long g
%     P = 100 ; 
%     k = 49 ; 
%     A = randn(P,k) ; 
%     a = randn(P,1) ;
%     Aast = inv(A'*A) ; 
%     B = [A a] ; 
%     BastReal = inv(B'*B) ;
%     b = randn(P,1) ; 
%     xold = Aast*(A'*b) ; 
%     r = b - A*xold ; 
%     xREAL = BastReal*(B'*b) ;
load('tmp3.mat')
    
 
end
c = A'*a ; 
d = Aast*c ; 
s = a'*a-c'*d ; 
Bast = [Aast + (d*d')/s  -d/s; -d'/s  1/s] ; 
v = a'*r/s ; 
x = [(xold -d*v ); v] ;

 
