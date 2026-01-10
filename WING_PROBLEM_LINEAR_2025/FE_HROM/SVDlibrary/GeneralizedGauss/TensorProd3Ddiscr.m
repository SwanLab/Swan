function [MPOINTS Mglo x w]  =TensorProd3Ddiscr(M,xLIM)

if nargin == 0
    M = [2,2,2];
    xLIM = [-1 1 ; -1 1; -1 1] ;
elseif    nargin ==1
   xLIM = [-1 1 ; -1 1; -1 1] ;
    
end

%
if M(1) ~=M(2)
    error('Not valid discretization  for tensor product rules ')
end
doffXLIM = xLIM-[-1 1 ; -1 1; -1 1];
if norm(doffXLIM)~=0
    error('Not valid limits (it should be a biunit square)')
end


[x1d, w1d] = GaussQuad(M(1), -1, 1) ;
%
%xxx =[]  ;
%yyy = [] ;
iacum = 1;
w = [] ;
for i = 1:length(x1d)
    for j = 1:length(w1d)
        for k = 1:length(w1d)
        %       xxx(iacum)  = x1d(i) ;
        %      yyy(iacum) =  x1d(j) ;
        w(iacum) = w1d(i)*w1d(j)*w1d(k) ;
        iacum = iacum+1 ;
        end
    end
end

x{1} = x1d' ;
x{2} = x1d' ;
x{3} = x1d' ; 
MPOINTS = M ;
Mglo = MPOINTS(1)*MPOINTS(2)*MPOINTS(3) ;
w = w' ;
