function [S,U]= PRANGLES(A,B,M,epsilon_min,epsilon_max)
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/DOMAINdecompositionGEN/RevampTheory.mlx
% JAHO, 11-Apr-2019 . Source: Golub + (My contribution: incorporation of the M matrix )
if nargin == 0
    load('tmp1.mat')
    A = B ;
    M = speye(size(M));
end

if ~isempty(M) 
    Mchol = chol(M) ;
    A = Mchol*A ;
    B = Mchol*B ;    
end

DATALOC.RELATIVE_SVD = 1 ;
TOL_SVD  = 1e-6 ;
[QA,SA,~] = SVDT(A,TOL_SVD,DATALOC) ;
[QB,SB,~] = SVDT(B,TOL_SVD,DATALOC) ;
[Y,S,Z] = SVDT(QA'*QB) ; 
r = find(((S > epsilon_min).*(S <= epsilon_max))==1);



if isempty(M)
    U = QA*Y(:,r) ;

else
    U = Mchol\(QA*Y(:,r)) ;

end

sA = size(QA,2) ;
sB = size(QB,2) ;

