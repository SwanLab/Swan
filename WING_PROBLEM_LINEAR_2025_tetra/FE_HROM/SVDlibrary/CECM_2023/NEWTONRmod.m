function [xNEW,wNEW,VARCnew,SALIRloc,SALIR,DATALOC] = NEWTONRmod(normB,VARCnew,wNEW,DATALOC,TOL,b,xNEW,VAR_SMOOTH_FE,POLYINFO)


SALIRloc = 0 ;
kMAX = DATALOC.MaxIterationsNR_ElimPoints  ;    k=1;    nF = 1e10 ;    
nF_old = nF ;    nITER_increase = 0 ;    maxITER_increase = DATALOC.maxITER_allowed_residual_withoutDECREASE; % = 4 ;
ISOUT = 0 ; SALIR = 0 ;
VARCnew.ListElementsInTransitionINNERloop = {} ;
while  nF>=TOL && k<= kMAX && nITER_increase<=maxITER_increase
    [xkp1, wkp1, nF, POLYINFO,ISOUT,VARCnew  ] = UPDATEPOSW(wNEW,b,xNEW,DATALOC,VAR_SMOOTH_FE,POLYINFO,VARCnew) ;
    nF = nF/normB ;    nd = (xkp1-xNEW).^2 ;    nd = max(sqrt(sum(nd,2))) ;
    disp(['Iteration k=',num2str(k),',  error residual =',num2str(nF),' MAX NORM incre DISPL =',num2str(nd),POLYINFO.MESSAGE_RANK]) ; 
    if ISOUT == 1
        break
    end
    DecreaseResidual = (nF_old-nF)/nF_old ;
    % Tolerance to consider that the residual is actually decreasing
    TOL_decrease_tolerance = 1e-6;
    if  DecreaseResidual < TOL_decrease_tolerance
        nITER_increase = nITER_increase + 1 ;
    end    
    xNEW = xkp1;
    wNEW = wkp1 ;
    k = k+1 ;
    nF_old = nF ;
end
DATALOC.TOL_low =  DATALOC.MaxIterationsNR_ElimPoints;
% ------------------------------------------------------
if ( ISOUT ==1  ) || (k>kMAX  &&  nF>DATALOC.TOL_low ) || nITER_increase >maxITER_increase
 %   iremove = iremove + 1;
 %   xNEW = xOLD ;
 %   wNEW = wOLD ;
    
  
    
  %  VARCnew = VARC ;
    SALIRloc = 1;
    
else
    SALIR = 1 ;
    
    % Points belonging to POINTSRp are moved to POINTSl
    VARCnew.POINTSl = [VARCnew.POINTSl(:); VARCnew.POINTSRp(:)] ;
    VARCnew.POINTSRp = [] ;
    
    
end
