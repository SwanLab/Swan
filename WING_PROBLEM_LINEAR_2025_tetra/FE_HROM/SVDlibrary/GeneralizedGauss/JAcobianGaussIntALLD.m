function D= JAcobianGaussIntALLD(dPHIk_y,xNEW,wNEW,PHIk_y,xINT,VSinv,DATALOC)

if nargin == 0
    load('tmp1.mat')
    
end
% Derivatives at points x(zNEW) ---quadratic interpolation
D = [] ;
% %dbstop('6')
% if    isempty(DATALOC.FUN_GRADIENT_DIRECTLY)
%     for idim = 1:length(xINT)
%         DerPHIt = zeros(size(PHIk_y,2) ,size(xNEW,1)) ;
%         %for i = 1:length(zNEW)
%         for j = 1:size(PHIk_y,2)
%             dPHIloc = reshape(PHI_der{idim}(:,j),size(xINT{1},1),[]);
%             DerPHIt(j,:) =  interp2(xINT{1},xINT{2},dPHIloc,xNEW(:,1),xNEW(:,2),'cubic');
%
%         end
%         DerPHIt = bsxfun(@times,DerPHIt',wNEW)' ;
%         D =[D DerPHIt] ;
%     end
% else
%  %   [dummy dPHIk_y]= EvaluateBasisFunctionAtX(xNEW,VSinv,DATALOC,1,0);
%
%       if DATALOC.APPROX_FUN__DERI.ACTIVE == 0
%       [dummy dPHIk_y]= EvaluateBasisFunctionAtX(xNEW,VSinv,DATALOC,1,0);
%     else
%         % Approximated method, fun. and derivatives
%         % ------------------------------------------- April-2020
%         disp('Approx. fun. and der. ...')
%         tic
%          [dummy dPHIk_y_old]= EvaluateBasisFunctionAtX(xNEW,VSinv,DATALOC,1,0);
%         [dummy dPHIk_y]= EvaluateBasisFunctionAtX_approx(xNEW, DATALOC.APPROX_FUN__DERI)  ;
%
%          EEE1 = norm(dPHIk_y{1}-dPHIk_y_old{1},'fro')./norm(dPHIk_y_old{1},'fro')*100 ;
%          EEE2 = norm(dPHIk_y{2}-dPHIk_y_old{2},'fro')./norm(dPHIk_y_old{2},'fro')*100 ;
%         disp(['Error approx Dx = ',num2str(EEE1),' %'])
%
%         disp(['Error approx Dy = ',num2str(EEE1),' %'])
%         toc
%         disp('...Done')
%     end
%
%


 