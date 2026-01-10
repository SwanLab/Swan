function [xNEW,wNEW,CONVERGENCE,DATALOC] = Iter_RemPnt2Dapprox_END(xNEW,wNEW,xINT,PHI,PHI_der,b,DATALOC)

%dbstop('5')
if nargin == 0
    load('tmp1.mat')
    
    %   load('tmpPR.mat')
    %    xNEW = x_ProducRule
    %    wNEW = w_ProducRule
    
    load('tmpMINE.mat','wGAUSS','xGAUSS')
    xNEW = xGAUSS ;
    wNEW = wGAUSS ;
    
end
% Evaluation of the basis functions at the points xNEW


% if isempty(DATALOC.FUN_GRADIENT_DIRECTLY)
%     PHIk_y = zeros(size(xNEW,1),size(PHI,2));
%
%     for i=1:size(PHI,2)
%         PHIloc = reshape(PHI(:,i),size(xINT{1},1),[]);
%         PHIk_y(:,i) = interp2(xINT{1},xINT{2},PHIloc,xNEW(:,1),xNEW(:,2),'cubic') ;
%     end
% else
% dbstop('20')
%     if DATALOC.APPROX_FUN__DERI.ACTIVE == 0
%         PHIk_y = EvaluateBasisFunctionAtX(xNEW,VSinv,DATALOC,0,1)  ;
%     else
% Approximated method, fun. and derivatives
% ------------------------------------------- April-2020
% disp('Approx. fun. and der. ...')
%  tic
%  PHIk_y_old = EvaluateBasisFunctionAtX(xNEW,VSinv,DATALOC,0,1)  ;


PHIk_y = EvaluateBasisFunctionAtX_approx(xNEW, DATALOC.APPROX_FUN__DERI)  ;


%  EEE = norm(PHIk_y-PHIk_y_old,'fro')./norm(PHIk_y_old,'fro')*100 ;
%  disp(['Error approx = ',num2str(EEE),' %'])
%  toc
disp('...Done')
%   end
%end

errCOMP = norm(PHIk_y'*wNEW-b)./norm(b)    % This error measures how approximate is the SVD spatial decomp.

% ------------------T
% Criterion for selecting the point to be removed: significance index
% This process is repeated
disp('*********************************************************************************************')
disp([' CUB. RULE  with    ', num2str(length(xNEW)-1),' points ( of ',num2str(DATALOC.mINI),')'])
disp('**********************************************************************************************')
TOL = DATALOC.TOL ;
% Criterion for removing points
S_crit = CriterRemovePnts(DATALOC,wNEW,PHIk_y,b,PHI_der,xNEW,xINT);

[MMM indREM ]= sort(S_crit) ;

xOLD = xNEW ;
wOLD = wNEW ;
% ----------------------------------------------------
% New set of points (and weights)
m = size(xNEW,1) ;
p = size(PHI,2) ;
% Generalized variable
iremove = 1 ;
SALIR = 0 ;
%dbstop('43')
normB = norm(b);
%dbstop('37')
while SALIR == 0 & iremove <=size(xOLD,1)
    %----------------------------------------------------
    % Remove point indREM(iremove)
    xNEW(indREM(iremove),:) = [] ;
    wNEW(indREM(iremove)) = [] ;
    disp('----------------------------------------')
    disp(['Removing point =',num2str(indREM(iremove)),' (irem =',num2str(iremove),')'])
    disp('---------------------------------------')
    kMAX = DATALOC.kMAX  ;
    k=1;
    nF = 1e10 ;
    CONVERGENCE = 1;
    nF_old = nF ;
    nITER_increase = 0 ;
    maxITER_increase = 4 ;
    xORG = xNEW ;
    % wORG = wNEW ;
    
    while  nF>=TOL & k<= kMAX & nITER_increase<=maxITER_increase
        
        
        % dbstop('74')
        [xkp1 wkp1 nF ISNEGATIVE  ] = IterNR_GAuss2Dapprox_END(wNEW,PHI,b,xNEW,PHI_der,xINT,DATALOC) ;
        nF = nF/normB ;
        
        dX = abs(xkp1-xORG) ;
        
        dxREL = max(dX(:,1))/DATALOC.hx*100 ;
        dyREL = max(dX(:,2))/DATALOC.hy*100 ;
        
        
        maxDD = max(DATALOC.dxRELmax , DATALOC.dyRELmax )  ;
        disp(['Iteration k=',num2str(k),',  error =',num2str(nF),', dxMAX(%)=',num2str(dxREL), ', dyMAX(%)=',num2str(dyREL),...
            ' MAXdxy(%) =',num2str(maxDD)])
        
        % dbstop('73')
        %  if DATALOC.NEGATIVE_check_during_iterations==1
        if ISNEGATIVE ==1
            break
        end
        % else
        
        % end
        
        if nF>nF_old
            nITER_increase = nITER_increase + 1 ;
        end
        %%% NEW ITERATIOn
        %  disp('Iteration k=2')
        xNEW = xkp1;
        wNEW = wkp1 ;
        
        k = k+1 ;
        nF_old = nF ;
    end
    % ------------------------------------------------------
    
    if ISNEGATIVE ==1 | (k>kMAX  &  nF>DATALOC.TOL_low ) | nITER_increase >maxITER_increase
        % -----------
        %     dbstop('58')
        
        iremove = iremove + 1;
        xNEW = xOLD ;
        wNEW = wOLD ;
        
    else
        SALIR = 1 ;
        %    if k>kMAX
        %       CONVERGENCE = 0 ;
        
        
        %  end
        DATALOC.dxRELmax = max(DATALOC.dxRELmax,dxREL) ;
        DATALOC.dyRELmax = max(DATALOC.dyRELmax,dyREL) ;
        
        
        
    end
end

if  iremove > length(xOLD)
    CONVERGENCE = 0 ;
end
