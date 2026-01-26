function [R,ANGLES,sa,sb,coeff_A ]= IntersectionSubspaces(A,B,TOL_SVD,TOL_ANGLE)
% Intersection between the column space of two matrices A and B
% OUTPUT: Orthogonal matrix for the column space of the intersection 
% JAHO, 11-Apr-2019 . Source: Golub 
 

DATALOC.RELATIVE_SVD = 1 ;
[QA,SA,VA] = SVDT(A,TOL_SVD,DATALOC) ;
coeff_A = bsxfun(@times,VA',1./SA)' ;
[QB,SB,VB] = SVDT(B,TOL_SVD,DATALOC) ;
[Y,CT,Z] = SVDT(QA'*QB) ;
ANGLES = real(acos(CT))*180/pi ;
nREACT = length(find(ANGLES <= TOL_ANGLE));
R = QA*Y(:,1:nREACT) ;

% 
coeff_A = coeff_A*Y(:,1:nREACT) ;  % R = A*coeff_A

sa = size(QA,2); 
sb = size(QB,2); 