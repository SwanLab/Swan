function [xNEW,wNEW,CONVERGENCE,DATALOC,POLYINFO,VARCnew ] = IterativeSearchRESZERO(xNEW,wNEW,b,DATALOC,VAR_SMOOTH_FE,VARC)
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/ContinuousEmpiricalCubatureM/Paper_hernandez2021ecm/07_BEAMS_NOGRAD/ControlledVariables.mlx
% Adaptation of /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/LARGE_STRAINS/ContinuousECM/ControlPointsAlgLARGE.m
% so that it can handle removal of a given point in a stepwise fashion
% -------------------------------------------------------------------------------------
if nargin == 0
    load('tmp1.mat')
    DATALOC.SECOND_STAGE_ITERATIONS  =10;
end


DATALOC = DefaultField(DATALOC,'MSGPRINT',{}) ;

POLYINFO = [] ; % Information about local searchs (for FE_interpolation option)
%DATAINeval = [] ;
[PHIk_y,dPHIk_y,POLYINFO]=     EvaluateBasisFunctionALL(xNEW,DATALOC,VAR_SMOOTH_FE,POLYINFO)  ;

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

%ONLYweights
switch DATALOC.criterion_for_removing_weights
    case 'STANDARD'
        S_crit = wNEW(POINTS_F).*sum(PHIk_y(POINTS_F,:).*PHIk_y(POINTS_F,:),2) ;
    case 'ONLYweights'
        error('Use the standard option; it has proved to be  more reliable')
        % S_crit = wNEW;
end
[S_crit_sorted, indREM ]= sort(S_crit) ;

xOLD = xNEW ;
wOLD = wNEW ;
VARCnew = VARC ;


% Generalized variable
iremove = 1 ;
SALIR = 0 ;
%dbstop('43')
normB = norm(b);
CONVERGENCE = 1;


%
while SALIR == 0 && iremove <=length(S_crit)
    iterWEIGHTS = 1 ;
    wcontr_old = 0 ;
    ipoint_control = 0 ;  SALIRloc = 0 ;
    
    % Me move the indexes of points with fixed positions to the set of
    % indexes with points free to change positions and weights
    
    
    if isempty(VARC.POINTSRpFIXED)
        VARCnew.POINTSl = [VARCnew.POINTSl(:); VARCnew.POINTSRp(:)] ;
        VARCnew.POINTSRp = [] ;
        
    else
        VARCnew_POINTSRp_old = VARCnew.POINTSRp ;
        [VARCnew.POINTSRp,~,IndexesP] =   intersect(VARCnew.POINTSRpFIXED,VARCnew_POINTSRp_old) ;
        % Points which are not in the intersection are moved to VARCnew.POINTSl(:)
        IndexMoveToL = setdiff(1:length(VARCnew_POINTSRp_old),IndexesP) ;
        VARCnew.POINTSl = [VARCnew.POINTSl(:); VARCnew_POINTSRp_old(IndexMoveToL)] ;
    end
    VARCnew.ListFirstElementJumpIteration = {} ;
    
     DATALOC.HISTORY_LOCAL.x = {} ;
     DATALOC.HISTORY_LOCAL.w  = {} ;
     DATALOC.HISTORY_LOCAL.ControlPoints = [] ;
    
    while iterWEIGHTS <=  DATALOC.SECOND_STAGE_ITERATIONS && SALIRloc == 0
        % This is the loop in which the weight of the "control" point is
        % gradually diminished
        [xNEW,wNEW,VARCnew,wcontr_old,ipoint_control,nF,ISNEGATIVE,ISOUT,SALIRloc,SALIR,iremove,DATALOC ] =...
            NewtonSearchLoc(normB,indREM,iremove,VARCnew,wNEW,DATALOC,TOL,b,xNEW,VAR_SMOOTH_FE,POLYINFO,iterWEIGHTS,...
            wcontr_old,ipoint_control,...
            xOLD,wOLD,VARC) ;
       
        if SALIRloc == 0
            DATALOC.HISTORY_LOCAL.x{end+1} = xNEW ;
            DATALOC.HISTORY_LOCAL.w{end+1}  = wNEW ;
            if iterWEIGHTS == 1
            DATALOC.HISTORY_LOCAL.ControlPoints= ipoint_control ; 
            end
        end
         iterWEIGHTS = iterWEIGHTS + 1;
    end
    
end
if  iremove >length(S_crit)
    CONVERGENCE = 0 ;
else
    DATALOC.REMOVED_INDEX = indREM(iremove) ;
end
