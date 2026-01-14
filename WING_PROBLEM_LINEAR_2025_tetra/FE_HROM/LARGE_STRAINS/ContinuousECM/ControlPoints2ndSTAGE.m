function  [xOLD,wOLD,DATALOC,POLYINFO,VARC,HISTORY,ELEMENTS_TO_PLOT] = ControlPoints2ndSTAGE(xOLD,wOLD,b,DATALOC,VAR_SMOOTH_FE,VARC,HISTORY,ELEMENTS_TO_PLOT)
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/ContinuousEmpiricalCubatureM/Paper_hernandez2021ecm/07_BEAMS_NOGRAD/ControlledVariables.mlx
if nargin == 0
    load('tmp_216.mat')
    %  DATALOC.SECOND_STAGE_ITERATIONS = 20;
    %     DATALOC.EXCLUDE_ELEMENT_TRANSITION_LIST =[] ;
    %  DATALOC.EXCLUDE_ELEMENTS_TRIGGERED_NONCONVERGENCE_ITERATION = 1 ; %2;
end

if ~isempty(DATALOC.EXCLUDE_ELEMENT_TRANSITION_LIST) && DATALOC.Method_Evaluation_Basis_Integrand ~=2
    if  ischar(DATALOC.EXCLUDE_ELEMENT_TRANSITION_LIST)
        ListElements = load(DATALOC.EXCLUDE_ELEMENT_TRANSITION_LIST) ;
        DATALOC.EXCLUDE_ELEMENT_TRANSITION_LIST = ListElements(:,1) ;
    end
end
DATALOC.FORCE_POINTS_TO_REMAIN_INSIDE =1 ; % We use always this option here !

CONVERGENCE = 1 ;
iter = 1;

if isempty(VARC.POINTSRpFIXED)
    VARC.POINTSl = find(wOLD ~=0);% These are the indexes of the points for which both weights and positions are unknowns
    % (candidates for being eliminated)
    VARC.POINTSRp = []; % These are the indexes of the points for which   positions   are constrained, but weights aren't
    VARC.POINTSRpw = find(wOLD ==0); % These are the indexes of the points for which   weights and positions are constrained
else
    % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/ContinuousEmpiricalCubatureM/Paper_hernandez2021ecm/07_BEAMS_NOGRAD/ControlledVariables3D.mlx
    
    VARC_POINTSRp_old = VARC.POINTSRp ;
    [VARC.POINTSRp,~,IndexesP] =   intersect(VARC.POINTSRpFIXED,VARC_POINTSRp_old) ;
    % Points which are not in the intersection are moved to VARCnew.POINTSl(:)
    IndexMoveToL = setdiff(1:length(VARC_POINTSRp_old),IndexesP) ;
    VARC.POINTSl = [VARC.POINTSl(:); VARC_POINTSRp_old(IndexMoveToL)] ;
    
    
end

DATALOC.BenignJumps = [] ;




DATALOC = DefaultField(DATALOC,'Include2ndStageIterations_PlotEvolutionWeights',0) ; % 15-Oct-2022

HISTORY.INDEXES_FIRST_STAGE = 1:length(HISTORY.POINTS_all) ;
HISTORY.INDEXEX_SECOND_STAGE = {};
HISTORY.CONTROL_POINTS = [] ; 
%HISTORY.CONTROL_POINTS =[] ;
NumberIterFirstStage = length(HISTORY.INDEXES_FIRST_STAGE) ;
CurrentITERglo = NumberIterFirstStage  ;
% ------------------------------------------------------
while CONVERGENCE ==1  %&& length(xOLD)>1
    DATALOC.iter = iter ;
    %  [xNEW,wNEW,CONVERGENCE,DATALOC,POLYINFO,VARC] = ControlPointsAlgLARGE(xOLD,wOLD,b,DATALOC,VAR_SMOOTH_FE,VARC) ;
    [xNEW,wNEW,CONVERGENCE,DATALOC,POLYINFO,VARC ] = IterativeSearchRESZERO(xOLD,wOLD,b,DATALOC,VAR_SMOOTH_FE,VARC) ;
    if CONVERGENCE == 1
        
        [ELEMENTS_TO_PLOT,HISTORY,INDnneg,CurrentITERglo] = HistoryPointsUpdate(wNEW,xNEW,ELEMENTS_TO_PLOT,HISTORY,DATALOC,POLYINFO,iter,CurrentITERglo) ;
        DATALOC.HISTORY_LOCAL.x = {} ;
        DATALOC.HISTORY_LOCAL.w  = {} ;
        iter = iter + 1 ;
        xOLD = xNEW ;        wOLD = wNEW ;
    end
end
%HISTORY.CONTROL_POINTS = DATALOC.HISTORY_LOCAL.ControlPoints ;

if iter == 1
    INDnneg = find(wOLD ~=0);
end

disp('------------------------------------------------------------------------------------------------')
disp(['SECOND STAGE: Integration rule with ',num2str(length(INDnneg)),' of ',num2str(length(wNEW)),' points'])
disp('------------------------------------------------------------------------------------------------')
disp(['List of elements whose points are constraint to move across each other'])
disp(num2str(DATALOC.EXCLUDE_ELEMENT_TRANSITION_LIST'))
