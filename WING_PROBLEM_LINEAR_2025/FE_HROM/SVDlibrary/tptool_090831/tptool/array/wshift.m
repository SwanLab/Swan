function y = wshift(type,x,p)
%WSHIFT Shift Vector or Matrix.
%   Y = WSHIFT(TYPE,X,P) with TYPE = {1,'1','1d' or '1D'}
%   performs a P-circular shift of vector X.
%   The shift P must be an integer, positive for right to left
%   shift and negative for left to right shift.
%
%   Y = WSHIFT(TYPE,X,P) with TYPE = {2,'2','2d' or '2D'}
%   performs a P-circular shift of matrix X.
%   The shifts P must be integers. P(1) is the shift for rows
%   and P(2) is the shift for columns.
%
%   WSHIFT('1D',X) is equivalent to WSHIFT('1D',X,1)
%   WSHIFT('2D',X) is equivalent to WSHIFT('2D',X,[1 1])
%
%   Example 1D:
%     x = [1 2 3 4 5]
%     wshift('1D',x,1)  = [2 3 4 5 1]
%     wshift('1D',x,-1) = [5 1 2 3 4]
%
%   Example 2D:
%     x = [1 2 3;5 6 7]
%     wshift('2D',x,[1 1])  = [6 7 5;2 3 1]
%     wshift('2D',x,[-1,0]) = [5 6 7;1 2 3]

%   M. Misiti, Y. Misiti, G. Oppenheim, J.M. Poggi 01-Dec-97.
%   Last Revision: 01-May-1998.
%   Copyright (c) 1995-98 by The MathWorks, Inc.
%   $Revision: 1.2 $  $Date: 1998/07/16 14:33:59 $

if isempty(x) | all(p==0) , y = x; return; end
switch type
  case {1,'1','1d','1D'}
    if nargin<3 , p = 1; end
    L = length(x);
    p = rem(p,L);
    if p<0 , p = L+p; end
    y = x([p+1:L,1:p]);

  case {2,'2','2d','2D'}
    if nargin<3 , p = [1 1]; end
    S = size(x);
    p = rem(p,S);
    k = (p<0);
    p(k) = S(k)+p(k);
    y = x([p(1)+1:S(1),1:p(1)],[p(2)+1:S(2),1:p(2)]);
end
