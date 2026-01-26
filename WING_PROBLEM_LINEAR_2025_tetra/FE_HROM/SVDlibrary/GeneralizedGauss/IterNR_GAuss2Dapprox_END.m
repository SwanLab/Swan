function [xNEW wNEW nF ISNEGATIVE  ] = IterNR_GAuss2Dapprox_END(wNEW,PHI,b,xNEW,PHI_der,xINT,DATALOC)

%dbstop('4')
if nargin == 0
    load('tmp2.mat')
end

ISNEGATIVE = 0 ; 

% if isempty(DATALOC.FUN_GRADIENT_DIRECTLY)
%     PHIk_y = zeros(length(wNEW),size(PHI,2));
%     
%     for i=1:size(PHI,2)
%         PHIloc = reshape(PHI(:,i),size(xINT{1}IterNR_GAuss2Dapprox,1),[]);
%         PHIk_y(:,i) = interp2(xINT{1},xINT{2},PHIloc,xNEW(:,1),xNEW(:,2),'cubic') ;
%     end
% else
 %   PHIk_y = EvaluateBasisFunctionAtX(xNEW,VSinv,DATALOC,0,1)  ;
    
%      if DATALOC.APPROX_FUN__DERI.ACTIVE == 0
%         PHIk_y = EvaluateBasisFunctionAtX(xNEW,VSinv,DATALOC,0,1)  ;
%     else
        % Approximated method, fun. and derivatives
        % ------------------------------------------- April-2020
    %    disp('Approx. fun. and der. ...')
     %   tic
      %   PHIk_y_old = EvaluateBasisFunctionAtX(xNEW,VSinv,DATALOC,0,1)  ;
        PHIk_y = EvaluateBasisFunctionAtX_approx(xNEW, DATALOC.APPROX_FUN__DERI)  ;
       % EEE = norm(PHIk_y-PHIk_y_old,'fro')./norm(PHIk_y_old,'fro')*100 ; 
        %disp(['Error approx = ',num2str(EEE),' %'])
        %toc
        %disp('...Done')
%     end
%     
%     
% end
 
m = length(wNEW);
bNEW = PHIk_y'*wNEW ; 
% Equation to be solved
% F =  (b-bNEW)
% Residual

%disp(' Residual norm(Fk)')
Fk = b-bNEW ; 
nF = norm(Fk)  ;

% Computation of the Jacobian matrix 
D= JAcobGaussInt2Dapprox(PHI_der,xNEW,wNEW,PHIk_y,xINT,DATALOC) ;

% Updated direction
delta_q = D\Fk ;

% dxMAX = max(abs(delta_q(1:m))) ;
% dyMAX =  max(abs(delta_q(m+1:2*m))) ;


%New point 
 q_kp1 = [xNEW(:);wNEW] +delta_q ;
 
 
 
wNEW = q_kp1(2*m+1:end) ; 
xNEW = [q_kp1(1:m) , q_kp1(m+1:2*m) ] ; 

 
xxLIM =DATALOC.xLIM(1,:) ;
yyLIM = DATALOC.xLIM(2,:) ;


% xxLIM = [min(min(xINT{1}))  max(max(xINT{1})) ] ; 
% yyLIM = [min(min(xINT{2}))  max(max(xINT{2})) ] ; 
 

% if all(wNEW>=0) & all(xNEW(:,1) >= xxLIM(1)) & all(xNEW(:,1) <= xxLIM(2)) & ...
%         all(xNEW(:,2) >= yyLIM(1)) & all(xNEW(:,2) <= yyLIM(2))
%    % disp('All points are admissible !!!!!!!')
% else
%   %dbstop('34')
%     disp('Inadmissible  point ...')
%     ISNEGATIVE = 1 ; 
% end

