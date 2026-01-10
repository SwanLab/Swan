function D = eig3(A)
% function D = eig3(A)
%
% Compute in one shot the eigen-values of multiples (3 x 3) matrices
%
% INPUT:
%   A: (3 x 3 x n) array
% OUTPUT:
%   D: (3 x n). EIG3 returns in D(:,k) three eigen-values of A(:,:,k)
%
% See also: CardanRoots, eig2, eig
%
% Author: Bruno Luong <brunoluong@yahoo.com>
% History:
%     Original 20-May-2010

if size(A,1) ~= 3 || size(A,2) ~= 3
    error('A must be [3x3xn] array');
end

A = reshape(A, 9, []).';

P3 = 1;
% Trace
P2 = -(A(:,1)+A(:,5)+A(:,9));
% Principal minors
M11 = A(:,5).*A(:,9) - A(:,8).*A(:,6);
M22 = A(:,9).*A(:,1) - A(:,3).*A(:,7);
M33 = A(:,1).*A(:,5) - A(:,4).*A(:,2);
P1 = (M11 + M22 + M33);
% Determinant
P0 = - A(:,1).*M11 ...
    + A(:,4).*(A(:,2).*A(:,9)-A(:,8).*A(:,3)) ...
    - A(:,7).*(A(:,2).*A(:,6)-A(:,5).*A(:,3));

D = CardanRoots(P3, P2, P1, P0).';

end
