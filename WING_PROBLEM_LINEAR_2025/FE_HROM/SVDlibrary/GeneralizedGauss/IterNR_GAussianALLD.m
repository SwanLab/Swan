function [xNEW wNEW nF ISNEGATIVE  ] = IterNR_GAussianALLD(wNEW,PHI,b,xNEW,PHI_der,xINT,VSinv,DATALOC)

%dbstop('4')
if nargin == 0
    load('tmp2.mat')
end

DATALOC = DefaultField(DATALOC,'USE_SVD_SOLVE_EQUATIONS',1) ;

 [PHIk_y,PHI_der] = EvaluateBasisFunctionAnalytical(xNEW,VSinv,DATALOC,1,1)  ;
    
m = length(wNEW);
bNEW = PHIk_y'*wNEW ;
% Equation to be solved
% F =  (b-bNEW)
% Residual
Fk = b-bNEW ;
nF = norm(Fk)  ;
% Computation of the Jacobian matrix
ndim = length(PHI_der) ; 
D = cell(1,ndim+1) ; 
for idim = 1:ndim
    D{idim} =  bsxfun(@times,PHI_der{idim},wNEW)' ; 
end
D{end} =  PHIk_y' ; 
D = cell2mat(D) ; 


 
% Updated direction
m = length(wNEW)  ; 
delta_q = SolutionUnderdeterminedSystem(DATALOC,D,Fk,wNEW,m) ;

 
if DATALOC.TOLERANCE_FORCE_POINTS_TO_REMAIN_WITHIN_THE_DOMAIN ==0
    % Standard option. If a point is outside the domain, the solution is
    % rejected (ISNEGATIVE ==1)
    [xNEW,wNEW,ISNEGATIVE] = UpdateSolutionWeightsPosition(xNEW,wNEW,delta_q,DATALOC)  ;
else 
    % Attempt to improve the performance of the method by forcing all
    % points to remain within the domain (Nov-2021)
    
    [xNEW,wNEW,ISNEGATIVE] = UpdateSolutionWeightsPositionFORCE(xNEW,wNEW,delta_q,DATALOC)  ;
end


% xxLIM = [min(min(xINT{1}))  max(max(xINT{1})) ] ;
% yyLIM = [min(min(xINT{2}))  max(max(xINT{2})) ] ;


% BEFORE OCTOBER-2021

% if all(wNEW>=0) & all(xNEW(:,1) >= xxLIM(1)) & all(xNEW(:,1) <= xxLIM(2)) & ...
%         all(xNEW(:,2) >= yyLIM(1)) & all(xNEW(:,2) <= yyLIM(2))
%    % disp('All points are admissible !!!!!!!')
% else
%   %dbstop('34')
%     disp('Inadmissible  point ...')
%     ISNEGATIVE = 1 ;
% end





