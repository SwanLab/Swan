function [xNEW wNEW nF ISNEGATIVE  ] = IterNR_GAussian2D(wNEW,PHI,b,xNEW,PHI_der,xINT,VSinv,DATALOC)

%dbstop('4')
if nargin == 0
    load('tmp2.mat')
end

DATALOC = DefaultField(DATALOC,'USE_SVD_SOLVE_EQUATIONS',1) ;

if isempty(DATALOC.FUN_GRADIENT_DIRECTLY)
    PHIk_y = zeros(length(wNEW),size(PHI,2));
    
    for i=1:size(PHI,2)
        PHIloc = reshape(PHI(:,i),size(xINT{1},1),[]);
        PHIk_y(:,i) = interp2(xINT{1},xINT{2},PHIloc,xNEW(:,1),xNEW(:,2),'cubic') ;
    end
else
    %   PHIk_y = EvaluateBasisFunctionAtX(xNEW,VSinv,DATALOC,0,1)  ;
    
    if DATALOC.APPROX_FUN__DERI.ACTIVE == 0
        PHIk_y = EvaluateBasisFunctionAtX(xNEW,VSinv,DATALOC,0,1)  ;
    else
        % Approximated method, fun. and derivatives
        % ------------------------------------------- April-2020
        disp('Approx. fun. and der. ...')
        tic
        PHIk_y_old = EvaluateBasisFunctionAtX(xNEW,VSinv,DATALOC,0,1)  ;
        PHIk_y = EvaluateBasisFunctionAtX_approx(xNEW, DATALOC.APPROX_FUN__DERI)  ;
        EEE = norm(PHIk_y-PHIk_y_old,'fro')./norm(PHIk_y_old,'fro')*100 ;
        disp(['Error approx = ',num2str(EEE),' %'])
        toc
        disp('...Done')
    end
    
    
end

m = length(wNEW);
bNEW = PHIk_y'*wNEW ;
% Equation to be solved
% F =  (b-bNEW)
% Residual

%disp(' Residual norm(Fk)')
Fk = b-bNEW ;
nF = norm(Fk)  ;

% Computation of the Jacobian matrix
D= JAcobianGaussInt2D(PHI_der,xNEW,wNEW,PHIk_y,xINT,VSinv,DATALOC) ;

% Updated direction
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





