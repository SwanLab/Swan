function [R,ANGLES,sA,sB ]= IntersectionSubspaces_M(A,B,M,TOL_SVD,TOL_ANGLE)
% Intersection between the column space of two matrices A and B
% OUTPUT: M-Orthogonal matrix for the column space of the intersection 
% JAHO, 11-Apr-2019 . Source: Golub + (My contribution: incorporation of the M matrix ) 
if nargin == 0
    load('tmp1.mat')
    A = B ; 
    M = speye(size(M));
end
Mchol = chol(M) ;
A = Mchol*A ;
B = Mchol*B ;

DATALOC.RELATIVE_SVD = 1 ;
[QA,SA,~] = SVDT(A,TOL_SVD,DATALOC) ;
[QB,SB,~] = SVDT(B,TOL_SVD,DATALOC) ;
[Y,CT,Z] = SVDT(QA'*QB) ;
ANGLES = real(acos(CT))*180/pi ;
nREACT = length(find(ANGLES <= TOL_ANGLE));
R = Mchol\(QA*Y(:,1:nREACT)) ;

sA = size(QA,2) ; 
sB = size(QB,2) ; 

 