function [xNEW,wNEW,CONVERGENCE,DATALOC,POLYINFO,VARCnew ] = MAKE1ZERO(xOLD,wOLD,b,DATALOC,VAR_SMOOTH_FE,VARC,POLYINFO)
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/ContinuousEmpiricalCubatureM/Paper_hernandez2021ecm/07_BEAMS_NOGRAD/ControlledVariables.mlx
% Adaptation of /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/LARGE_STRAINS/ContinuousECM/ControlPointsAlgLARGE.m
% so that it can handle removal of a given point in a stepwise fashion
% -------------------------------------------------------------------------------------
if nargin == 0
    load('tmp.mat')
    
end

[U_X,~,POLYINFO]=     EVALBASIS(xOLD,DATALOC,VAR_SMOOTH_FE,POLYINFO)  ;
% ------------------T
% Criterion for selecting the point to be removed: significance index
% This process is repeated
POINTS_F = [VARC.POINTSl(:); VARC.POINTSRp(:)] ; % All candidate points

disp('*********************************************************************************************')
disp([' CUB. RULE  with    ', num2str(length(POINTS_F)-1),' points ( of ',num2str(DATALOC.mINI),')'])
disp('**********************************************************************************************')
TOL = DATALOC.TOL_NewtonRaphson_EliminationPoints ;
%-------------------------------
% Criterion for removing points
% ------------------------------
% WE can only move points within the set VARC.POINTS_L
S_crit = wOLD(POINTS_F).*sum(U_X(POINTS_F,:).*U_X(POINTS_F,:),2) ;
[~, indREM ]= sort(S_crit) ;


% Generalized variable
iremove = 1 ;
SALIR = 0 ;
%dbstop('43')
normB = norm(b);
CONVERGENCE = 1;



while SALIR == 0 && iremove <=length(S_crit)
    
    %----------------------------------------------------
    iremovLOC = indREM(iremove) ; % Indexes of point belonging to  POINTS_F which is constrained now (because we are going to set its weight to zero,
    
    
    %wcontr_old = 0 ;
    %ipoint_control = 0 ;
    %  SALIRloc = 0 ;
    % Me move the indexes of points with fixed positions to the set of
    % indexes with points free to change positions and weights
    VARC.POINTSl = [VARC.POINTSl(:); VARC.POINTSRp(:)] ;
    VARC.POINTSRp = [] ;
    DATALOC.HISTORY_LOCAL.x = {} ;
    DATALOC.HISTORY_LOCAL.w  = {} ;
    DATALOC.HISTORY_LOCAL.ControlPoints = [] ;
    
    [xNEW,wNEW,VARCnew,SALIRloc,SALIR,DATALOC ] =...
        SOLVERES(normB,iremovLOC,wOLD,DATALOC,TOL,b,xOLD,VAR_SMOOTH_FE,...
        POLYINFO,VARC,iremove)  ;
    
    
    
    if SALIRloc == 1
        iremove = iremove + 1;
    end
    
    
end
if  iremove >length(S_crit)
    CONVERGENCE = 0 ;
else
    DATALOC.REMOVED_INDEX = indREM(iremove) ;
end
