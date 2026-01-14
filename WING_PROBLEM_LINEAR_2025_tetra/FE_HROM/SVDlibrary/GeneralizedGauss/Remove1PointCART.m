function [xNEW,wNEW,CONVERGENCE,DATALOC] = Remove1PointCART(xNEW,wNEW,b,DATALOC,...
    VSinv,DATAFITTING)

%dbstop('5')
if nargin == 0
    load('tmp2.mat')
end

% Just for checking purpuses...
EVALUATE_GRADIENT = 0;
g = Evaluate_Basis_Grad_Analytic(xNEW,VSinv,DATALOC,EVALUATE_GRADIENT,DATAFITTING)  ;
errCOMP = norm(g*wNEW-b)./norm(b);

% ------------------
% Criterion for selecting the point to be removed: significance index
% This process is repeated
disp('*********************************************************************************************')
disp([' CUB. RULE  with    ', num2str(length(xNEW)-1),' points ( of ',num2str(DATALOC.mINI),')'])
disp('**********************************************************************************************')
DATALOC = DefaultField(DATALOC,'TOL_NewtonRaphson',1e-10) ;

TOL = DATALOC.TOL_NewtonRaphson ;
kMAX = 30 ;
DATALOC = DefaultField(DATALOC,'maxITER_increase_NewtonRaphson',10);
maxITER_increase = DATALOC.maxITER_increase_NewtonRaphson ;

% Criterion for removing points
%S_crit = wNEW.*sum(g.*g,1) ; % Before 29-Dec-2021 
S_crit = wNEW.*(sum(g.*g,1)') ;

[MMM indREM ]= sort(S_crit) ;

xOLD = xNEW ; wOLD = wNEW ; zOLD = DATALOC.zOLD ; zNEW = zOLD ;
% ----------------------------------------------------
iremove = 1 ; SALIR = 0 ; normB = norm(b);
while SALIR == 0 && iremove <=size(xOLD,1)
    %----------------------------------------------------
    xNEW(indREM(iremove),:) = [] ;     zNEW(indREM(iremove)) = [] ;     wNEW(indREM(iremove)) = [] ;
    disp('----------------------------------------')
    disp(['Removing point =',num2str(indREM(iremove)),' (irem =',num2str(iremove),')'])
    disp('---------------------------------------')
    k=1;     nF = 1e10 ;     CONVERGENCE = 1;
    nF_old = nF ;    nITER_increase = 0 ;
    
    while  nF>=TOL && k<= kMAX && nITER_increase<=maxITER_increase
        
        [xkp1, wkp1, nF, ISNEGATIVE  ] = UpdatePositionsWeights(wNEW,b,xNEW,VSinv,DATALOC,DATAFITTING) ;
        nF = nF/normB ;
        
        disp(['Iteration k=',num2str(k),',  error =',num2str(nF)])
        
       
        
        if ISNEGATIVE ==1
            break
        end
      %  TOL = 0.1*eps(nF_old);
         if nF>=nF_old % Change 5th-Nov-2021
%         % See Remove1PointCART_aux.mlx
      %   TOL_LOC = 1e-10 ;  ;
      %   dN = nF-nF_old ;
      %   if dN >0 || abs(dN)/(nF_old)<=TOL_LOC % Change 10th-Nov-2021
            nITER_increase = nITER_increase + 1 ;
       %  else
        %     nITER_increase = 0 ;  % Change 1-Dec-2021
        end
        xNEW = xkp1;
        wNEW = wkp1 ;
        k = k+1 ;
        nF_old = nF ;
    end
    % ------------------------------------------------------
    
    %if ISNEGATIVE ==1 || (k>kMAX  &&  nF>DATALOC.TOL_low ) ||
    %nITER_increase >maxITER_increase. Removed 10th-Nov-2021
    if ISNEGATIVE ==1 || k>kMAX   || nITER_increase >maxITER_increase
        iremove = iremove + 1;
        xNEW = xOLD ;
        wNEW = wOLD ;
        zNEW = zOLD ;
        
    else
        SALIR = 1 ;
        DATALOC.zNEW = zNEW ;
        
    end
end

if  iremove > length(xOLD)
    CONVERGENCE = 0 ;
end
