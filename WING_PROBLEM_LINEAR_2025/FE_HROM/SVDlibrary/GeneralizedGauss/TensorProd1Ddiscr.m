function [MPOINTS Mglo x w]  =TensorProd1Ddiscr(M,xLIM)

if nargin ==0
    load('tmp.mat')
end



[x1d, w1d] = GaussQuad(M(1), -1, 1) ;
%
%xxx =[]  ;
%yyy = [] ;
iacum = 1;
w = [] ;
for i = 1:length(x1d)
    for j = 1:length(w1d)
        %       xxx(iacum)  = x1d(i) ;
        %      yyy(iacum) =  x1d(j) ;
        w(iacum) = w1d(i)*w1d(j) ;
        iacum = iacum+1 ;
    end
end

x{1} = x1d' ;
x{2} = x1d' ;
MPOINTS = M ;
Mglo = MPOINTS(1)*MPOINTS(2) ;
w = w' ;
