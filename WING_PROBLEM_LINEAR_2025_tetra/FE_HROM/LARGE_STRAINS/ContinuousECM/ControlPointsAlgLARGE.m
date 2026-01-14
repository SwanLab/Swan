function [xNEW,wNEW,CONVERGENCE,DATALOC,POLYINFO,VARCnew] = ControlPointsAlgLARGE(xNEW,wNEW,b,DATALOC,VAR_SMOOTH_FE,VARC)
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/ContinuousEmpiricalCubatureM/Paper_hernandez2021ecm/07_BEAMS_NOGRAD/ControlledVariables.mlx
if nargin == 0
    load('tmp1.mat')
end


DATALOC = DefaultField(DATALOC,'MSGPRINT',{}) ;

POLYINFO = [] ; % Information about local searchs (for FE_interpolation option)
%DATAINeval = [] ;
[PHIk_y,dPHIk_y,POLYINFO]=     EvaluateBasisFunctionALL(xNEW,DATALOC,VAR_SMOOTH_FE,POLYINFO)  ;

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
POINTS_F = [VARC.POINTSl(:); VARC.POINTSRp(:)] ; % All candidate points

%  if length(POINTS_F)-1 ==47
%      disp('Borra esto')
%  end

disp('*********************************************************************************************')
disp([' CUB. RULE  with    ', num2str(length(POINTS_F)-1),' points ( of ',num2str(DATALOC.mINI),')'])
disp('**********************************************************************************************')
DATALOC = DefaultField(DATALOC,'TOL_NewtonRaphson_EliminationPoints',1e-8) ;
TOL = DATALOC.TOL_NewtonRaphson_EliminationPoints ;
%-------------------------------
% Criterion for removing points
% ------------------------------
% WE can only move points within the set VARC.POINTS_L

%ONLYweights
switch DATALOC.criterion_for_removing_weights
    case 'STANDARD'
        S_crit = wNEW(POINTS_F).*sum(PHIk_y(POINTS_F,:).*PHIk_y(POINTS_F,:),2) ;
    case 'ONLYweights'
        error('Use the standard option; it has proved to be  more reliable')
        % S_crit = wNEW;
end
[S_crit_sorted, indREM ]= sort(S_crit) ;
% weigSORTED = wNEW(POINTS_F(indREM));
% weigUNIFORM = sum(wNEW)/length(indREM) ;


%VARC.POINTS_F  = POINTS_F ;
xOLD = xNEW ;
wOLD = wNEW ;
VARCnew = VARC ;
% VARCnew.PHIk_y = PHIk_y ;
%  VARCnew.dPHIk_y = dPHIk_y ;


% ----------------------------------------------------
% Generalized variable
iremove = 1 ;
SALIR = 0 ;
%dbstop('43')
normB = norm(b);


%
while SALIR == 0 && iremove <=length(S_crit)
    
    % Elimination strategy
    % --------------------
    %     factorWEIGHT = weigSORTED(iremove)/weigUNIFORM ;
    %     disp(['factor WEIGHT =',num2str(factorWEIGHT)])
    %      disp(['factor Scrit =',num2str(S_crit_sorted(iremove))])
    
    if isempty(VARCnew.POINTSRpFIXED)
        % Approach before  11-Jan-2022
        [VARCnew,xNEW,wNEW ]= EliminationStrategy(indREM,iremove,VARCnew,xNEW,wNEW) ;
        
    else
        %
        [VARCnew,xNEW,wNEW ]= EliminationStrategyFIXED(indREM,iremove,VARCnew,xNEW,wNEW) ;
    end
    
    DATALOC = DefaultField(DATALOC,'MaxIterationsNR_ElimPoints',40) ;
    
    kMAX = DATALOC.MaxIterationsNR_ElimPoints  ;    k=1;    nF = 1e10 ;    CONVERGENCE = 1;
    nF_old = nF ;    nITER_increase = 0 ;    maxITER_increase = 4 ;
    ISOUT = 0 ;
    VARCnew.ListElementsInTransitionINNERloop = {} ;
    while  nF>=TOL & k<= kMAX & nITER_increase<=maxITER_increase
        
        [xkp1, wkp1, nF, ISNEGATIVE,POLYINFO,ISOUT,VARCnew  ] = ...
            UpdateCoordinatesPoints_CONTROL(wNEW,b,xNEW,DATALOC,VAR_SMOOTH_FE,POLYINFO,VARCnew) ;
        
        nF = nF/normB ;    nd = (xkp1-xNEW).^2 ;    nd = max(sqrt(sum(nd,2))) ;
        disp(['Iteration k=',num2str(k),',  error residual =',num2str(nF),' MAX NORM incre DISPL =',num2str(nd)])
        
        if ISOUT == 1
            break
        end
        %         if nF>nF_old % bEFORE 20-oCT-2022 % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/ContinuousEmpiricalCubatureM/Paper_hernandez2021ecm/09_ASSESSMENTpaper/Fun3D/README_Cubature3D.mlx
        %             nITER_increase = nITER_increase + 1 ;
        %         end
        
        DecreaseResidual = (nF_old-nF)/nF_old ;
        % Tolerance to consider that the residual is actually decreasing
        TOL_decrease_tolerance = 1e-6; % Set it to zero 
        if  DecreaseResidual < TOL_decrease_tolerance
            nITER_increase = nITER_increase + 1 ;
        end
        
        
        
        xNEW = xkp1;
        wNEW = wkp1 ;
        k = k+1 ;
        nF_old = nF ;
    end
    % ------------------------------------------------------
    
    if ( ISOUT ==1  ) | (k>kMAX  &  nF>DATALOC.TOL_NewtonRaphson_EliminationPoints ) | nITER_increase >maxITER_increase
        iremove = iremove + 1;
        xNEW = xOLD ;
        wNEW = wOLD ;
        VARCnew = VARC ;
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
    end
end
if  iremove >length(S_crit)
    CONVERGENCE = 0 ;
else
    DATALOC.REMOVED_INDEX = indREM(iremove) ;
end
