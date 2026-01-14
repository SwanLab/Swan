function [xNEW,wNEW,CONVERGENCE,DATALOC,POLYINFO] = RemovingPointsAlgLARGE(xNEW,wNEW,b,DATALOC,VAR_SMOOTH_FE)

%dbstop('5')


DATALOC = DefaultField(DATALOC,'MSGPRINT',{}) ;

POLYINFO = [] ; % Information about local searchs (for FE_interpolation option)
DATAINeval = [] ;
[PHIk_y,~,POLYINFO]=     EvaluateBasisFunctionAtX_LARGE(xNEW,DATAINeval,VAR_SMOOTH_FE,POLYINFO)  ;

%ELEMENTS_CONTAINING_xNEW = POLYINFO.ELEMENTS_CONTAINING_xNEW ; 

errCOMP = norm(PHIk_y'*wNEW-b)./norm(b);      % MEASURE OF  INTERPOLATION ERROR !!! 
if DATALOC.iter == 1
    DATALOC.errorFITTING = errCOMP ;
    disp('--------------------------------------------') ;
    disp(['Integration appr. error at k= 0 -->',num2str(errCOMP*100),' %', '(It measures the quality in evaluating the function via   fitting)']);
end

% ------------------T
% Criterion for selecting the point to be removed: significance index
% This process is repeated
disp('*********************************************************************************************')
disp([' CUB. RULE  with    ', num2str(length(xNEW)-1),' points ( of ',num2str(DATALOC.mINI),')'])
disp('**********************************************************************************************')
TOL = DATALOC.TOL ;
% Criterion for removing points

%ONLYweights
switch DATALOC.criterion_for_removing_weights
    case 'STANDARD'
    S_crit = wNEW.*sum(PHIk_y.*PHIk_y,2) ;
    case 'ONLYweights'
        error('Use the other option; is more reliable')
        S_crit = wNEW; 
end
[MMM indREM ]= sort(S_crit) ;
xOLD = xNEW ;
wOLD = wNEW ;
%ELEMENTS_CONTAINING_xOLD = ELEMENTS_CONTAINING_xNEW ; 
% ----------------------------------------------------
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
 %   ELEMENTS_CONTAINING_xNEW(indREM(iremove)) = [] ; 
    
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
    ISOUT = 0 ;
    % wORG = wNEW ;
    while  nF>=TOL & k<= kMAX & nITER_increase<=maxITER_increase
        
        [xkp1, wkp1, nF, ISNEGATIVE,POLYINFO,ISOUT  ] = UpdateCoordinatesPoints_LARGE(wNEW,b,xNEW,DATALOC,VAR_SMOOTH_FE,POLYINFO) ;
       
        
        
        nF = nF/normB ;
         nd = (xkp1-xNEW).^2 ;
         nd = max(sqrt(sum(nd,2))) ; 
        disp(['Iteration k=',num2str(k),',  error residual =',num2str(nF),' MAX NORM incre DISPL =',num2str(nd)])
       
        
        
%         if any(wkp1<0)
%             disp('Borrar esto')
%         end
        
        %     if ISNEGATIVE ==1    % Before 24-april-2020
        if ISOUT == 1   % || (ISNEGATIVE ==1  && nF<TOL)
            
            break
        end
        if nF>nF_old
            nITER_increase = nITER_increase + 1 ;
        end
        
        xNEW = xkp1;
        wNEW = wkp1 ;
        k = k+1 ;
        nF_old = nF ;
    end
    % ------------------------------------------------------
    
   % if (ISNEGATIVE ==1 || ISOUT ==1  ) | (k>kMAX  &  nF>DATALOC.TOL_low ) | nITER_increase >maxITER_increase
        
   if ( ISOUT ==1  ) | (k>kMAX  &  nF>DATALOC.TOL_low ) | nITER_increase >maxITER_increase
        iremove = iremove + 1;
        xNEW = xOLD ;
        wNEW = wOLD ;
    else
        SALIR = 1 ;
    end
end

if  iremove > size(xOLD,1)
    CONVERGENCE = 0 ;
    
else
    DATALOC.REMOVED_INDEX = indREM(iremove) ;
    
end
