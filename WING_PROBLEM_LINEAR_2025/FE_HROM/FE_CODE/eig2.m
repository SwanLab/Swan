function D = eig2(A)
% function D = eig2(A)
%
% Compute in one shot the eigen-values of multiples (2 x 2) matrices
%
% INPUT:
%   A: (2 x 2 x n) array
% OUTPUT:
%   D: (2 x n). EIG2 returns in D(:,k) three eigen-values of A(:,:,k)
%
% See also: ParabolaRoots, eig3, eig
%
% Author: Bruno Luong <brunoluong@yahoo.com>
% History:
%     Original 27-May-2010

if size(A,1) ~= 2 || size(A,2) ~= 2
    error('A must be [3x3xn] array');
end

A = reshape(A, 4, []).';

P3 = 1;
% Trace
P2 = -(A(:,1)+A(:,4));

% Determinant
P1 = A(:,1).*A(:,4) - A(:,2).*A(:,3);

% Find the roots of characteristic polynomials
D = ParabolaRoots(P3, P2, P1).';

end % eig2

