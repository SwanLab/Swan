function [xNEW,wNEW,VARCnew,wcontr_old,ipoint_control,nF,ISNEGATIVE,ISOUT,SALIRloc,SALIR,iremove,DATALOC] =...
    NewtonSearchLoc(normB,indREM,iremove,VARCnew,wNEW,DATALOC,TOL,b,xNEW,VAR_SMOOTH_FE,POLYINFO,...
    iterWEIGHTS,wcontr_old,ipoint_control,xOLD,wOLD,VARC)

% Elimination strategy
% --------------------
% Here is where we choose the "control" point (at the very first
% iteration); we also set the weight corresponding at the current iteration
% (linearly scaled)
if isempty(VARCnew.POINTSRpFIXED)
    [VARCnew,xNEW,wNEW,wcontr_old,ipoint_control ]= EliminationStrategySTEP(indREM,iremove,VARCnew,xNEW,wNEW,DATALOC,iterWEIGHTS,wcontr_old,ipoint_control) ;
else
    [VARCnew,xNEW,wNEW,wcontr_old,ipoint_control ]= ...
        ElimStratSTEPfixed(indREM,iremove,VARCnew,xNEW,wNEW,DATALOC,iterWEIGHTS,wcontr_old,ipoint_control) ;    
end


SALIRloc = 0 ;
kMAX = DATALOC.MaxIterationsNR_ElimPoints  ;    k=1;    nF = 1e10 ;    CONVERGENCE = 1;
nF_old = nF ;    nITER_increase = 0 ;    maxITER_increase = DATALOC.maxITER_allowed_residual_withoutDECREASE; % = 4 ;
ISOUT = 0 ; SALIR = 0 ;
VARCnew.ListElementsInTransitionINNERloop = {} ;

while  nF>=TOL & k<= kMAX & nITER_increase<=maxITER_increase
    
    [xkp1, wkp1, nF, ISNEGATIVE,POLYINFO,ISOUT,VARCnew  ] = UpdateCoordinatesPoints_CONTROL(wNEW,b,xNEW,DATALOC,VAR_SMOOTH_FE,POLYINFO,VARCnew) ;
    
    nF = nF/normB ;    nd = (xkp1-xNEW).^2 ;    nd = max(sqrt(sum(nd,2))) ;
    disp(['Iteration k=',num2str(k),',  error residual =',num2str(nF),' MAX NORM incre DISPL =',num2str(nd)])
    
    if ISOUT == 1
        break
    end
%     if nF>nF_old    % bEFORE 20-oCT-2022 % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/ContinuousEmpiricalCubatureM/Paper_hernandez2021ecm/09_ASSESSMENTpaper/Fun3D/README_Cubature3D.mlx
%         nITER_increase = nITER_increase + 1 ;
%     end    
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
    iremove = iremove + 1;
    xNEW = xOLD ;
    wNEW = wOLD ;
    
    % CONVERGENCE HAS NOT BEEN ACHIEVED
    % LET US IDENTIFY, IF ANY, THE INTERELEMENT JUMPS THAT TRIGGERED
    % THE NON-CONVERGENCE
    DATALOC =  IdentifyCriticalElementsJump(DATALOC,iterWEIGHTS,VARCnew) ; % 16-Oct-2002 (modified, not revised)
    
    VARCnew = VARC ;
    SALIRloc = 1;
    
else
    SALIR = 1 ;
    if  isempty(VARCnew.POINTSRpFIXED)
        % Points belonging to POINTSRp are moved to POINTSl
        VARCnew.POINTSl = [VARCnew.POINTSl(:); VARCnew.POINTSRp(:)] ;
        VARCnew.POINTSRp = [] ;
    else
        % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/ContinuousEmpiricalCubatureM/Paper_hernandez2021ecm/07_BEAMS_NOGRAD/ControlledVariables3D.mlx
        % Intersection between VARCnew.POINTSRp and
        % VARCnew.POINTSRpFIXED
        VARCnew_POINTSRp_old = VARCnew.POINTSRp ;
        [VARCnew.POINTSRp,~,IndexesP] =   intersect(VARCnew.POINTSRpFIXED,VARCnew_POINTSRp_old) ;
        % Points which are not in the intersection are moved to VARCnew.POINTSl(:)
        IndexMoveToL = setdiff(1:length(VARCnew_POINTSRp_old),IndexesP) ;
        VARCnew.POINTSl = [VARCnew.POINTSl(:); VARCnew_POINTSRp_old(IndexMoveToL)] ;
        
    end
    DATALOC = UpdateBenignElements(VARCnew,DATALOC) ;
    
end
