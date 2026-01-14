function [xNEW, wNEW, nR, ISNEGATIVE  ] = UpdatePositionsWeights(wNEW,b,xNEW,VSinv,DATALOC,DATAFITTING)

%dbstop('4')
if nargin == 0
    load('tmp2.mat')
end

DATALOC = DefaultField(DATALOC,'USE_SVD_SOLVE_EQUATIONS',1) ;

EVALUATE_GRADIENT = 1;  
 [g,Grad_g ]= Evaluate_Basis_Grad_Analytic(xNEW,VSinv,DATALOC,EVALUATE_GRADIENT,DATAFITTING)  ;
 
bNEW = g*wNEW ;
% Equation to be solved
% F =  (b-bNEW)
% Residual
r = b-bNEW ;
nR = norm(r)  ;
% Computation of the Jacobian matrix
ndim = length(Grad_g) ; 
D = cell(1,ndim+1) ; 
for idim = 1:ndim
    D{idim} =  bsxfun(@times,Grad_g{idim}',wNEW)' ; 
end
D{end} =  g ; 
D = cell2mat(D) ; 

% Truncated SVD of Jacobian matrix
DATALOCSVD.RELATIVE_SVD = 1;
TOL = 1e-10 ;
[UU,SS,VV] = SVDT(D,TOL,DATALOCSVD) ;

% Solution of the underdetermined system of equations (more unknowns than equation)
if length(SS) == size(D,1)
    % No need to correct rank
    % -------------------------------------------------------
    delta_q = D\r ;    % Matlab produces a sparse solution (minimize  the l1 norm of delta_q)
else
    disp(['Incomplete rank: number of equations = ',num2str(size(D,1)),'; rank = ',num2str(length(SS)),' (TOLrel =',num2str(TOL)  ,')']);
    UU = bsxfun(@times,UU',1./SS)'  ;
    delta_q = VV'\(UU'*r) ;   
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Return points that at outside the domain 
 [xNEW,wNEW,ISNEGATIVE] = UpdateSolutionWeightsPositionFORCE(xNEW,wNEW,delta_q,DATALOC)  ;
 

