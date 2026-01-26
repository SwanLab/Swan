function [xNEW,wNEW,CONVERGENCE,DATALOC] = Iter_RemovePointALLDcart(xNEW,wNEW,xINT,PHI,PHI_der,b,DATALOC,...
    VSinv)

%dbstop('5')
if nargin == 0
    load('tmp2.mat')
end

PHIk_y = EvaluateBasisFunctionAnalytical(xNEW,VSinv,DATALOC,0,1)  ;
errCOMP = norm(PHIk_y'*wNEW-b)./norm(b)

% ------------------
% Criterion for selecting the point to be removed: significance index
% This process is repeated
disp('*********************************************************************************************')
disp([' CUB. RULE  with    ', num2str(length(xNEW)-1),' points ( of ',num2str(DATALOC.mINI),')'])
disp('**********************************************************************************************')
TOL = DATALOC.TOL ;
% Criterion for removing points
S_crit = wNEW.*sum(PHIk_y.*PHIk_y,2) ;

[MMM indREM ]= sort(S_crit) ;

xOLD = xNEW ; wOLD = wNEW ; zOLD = DATALOC.zOLD ; zNEW = zOLD ;
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
    xNEW(indREM(iremove),:) = [] ;     zNEW(indREM(iremove)) = [] ;     wNEW(indREM(iremove)) = [] ;
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
    
    while  nF>=TOL & k<= kMAX & nITER_increase<=maxITER_increase
       
        [xkp1 wkp1 nF ISNEGATIVE  ] = IterNR_GAussianALLD(wNEW,PHI,b,xNEW,PHI_der,xINT,VSinv,DATALOC) ;
        nF = nF/normB ;
        
        dX = abs(xkp1-xORG) ;
        
        %         dxREL = max(dX(:,1))/DATALOC.hx*100 ;
        %         dyREL = max(dX(:,2))/DATALOC.hy*100 ;
        
        
        %     maxDD = max(DATALOC.dxRELmax , DATALOC.dyRELmax )  ;
        %         disp(['Iteration k=',num2str(k),',  error =',num2str(nF),', dxMAX(%)=',num2str(dxREL), ', dyMAX(%)=',num2str(dyREL),...
        %             ' MAXdxy(%) =',num2str(maxDD)])
        
        disp(['Iteration k=',num2str(k),',  error =',num2str(nF)])
        
        % dbstop('73')
        %  if DATALOC.NEGATIVE_check_during_iterations==1
        if ISNEGATIVE ==1
            
            break
        end
        % else
        
        % end
        
        %if nF>nF_old
        
        if nF>=nF_old % Change 5th-Nov-2021
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
        zNEW = zOLD ;
        
    else
        SALIR = 1 ;
        DATALOC.zNEW = zNEW ;
        %    if k>kMAX
        %       CONVERGENCE = 0 ;
        
        
        %  end
        %         DATALOC.dxRELmax = max(DATALOC.dxRELmax,dxREL) ;
        %         DATALOC.dyRELmax = max(DATALOC.dyRELmax,dyREL) ;
        
        
        
    end
end

if  iremove > length(xOLD)
    CONVERGENCE = 0 ;
end
