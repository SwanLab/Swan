function [xGAUSS,wGAUSS,DATA_ECM,VAR_SMOOTH_FE,POLYINFO,AUXVAR] ...
    = ConstrainedCECM_local(COORg,w,DATA_ECM,coorECM,VAR_SMOOTH_FE)
%dbstop('6')
if nargin ==0
    load('tmp1.mat')
    %  DATA_ECM.FORCE_POINTS_TO_REMAIN_INSIDE = 1 ;
    %  DATA_ECM.SECOND_STAGE_ITERATIONS = 10 ;
    %  DATA_ECM.ONE_STAGE_ITERATIONS = 0 ;
end

b=   DATA_ECM.ExactIntegral; % PHI'*W ;
VAR_SMOOTH_FE.ExactIntegral = b;

% Starting Cubature rule
% ----------------------
wOLD =w ;
xOLD = coorECM ;
PATHiniPOINTS = cell(size(xOLD,1),size(xOLD,1)) ;
DATA_ECM.mINI = size(xOLD,1) ;

for ipoint = 1:size(xOLD,1)
    PATHiniPOINTS{ipoint,1 } = xOLD(ipoint,:) ;
end
%INDEX_xNEW_respect_XOLD =  1:size(xOLD,1) ;
% At each iteration of the ensuing loop, a point of xOLD is removed
% We shall asign to the set of initial points the labels 1,2....npoints
% At each iteration, we have to ascertain which are the points that remains
% in the setPlotPathPointsGauss2D
% For plotting affecgted elements
%[~,~,POLYINFOloc]=     EvaluateBasisFunctionDIRECTFIT(xOLD,[],VAR_SMOOTH_FE,[]) ;
[~,~,POLYINFOloc]=     EvaluateBasisFunctionALL(xOLD,DATA_ECM,VAR_SMOOTH_FE,[]) ;



ELEMENTS_TO_PLOT =POLYINFOloc.ELEMENTS_CONTAINING_xNEW  ;


%NEW variables --> 24-Apr-2020 -->
HISTORY.POINTS{1} = xOLD ;
HISTORY.WEIGHTS{1} = wOLD ;
HISTORY.POINTS_all{1} = xOLD ;
HISTORY.WEIGHTS_all{1} = wOLD ;
HISTORY.ISALLPOSITIVE = 1  ;
HISTORY.ELEMENTS_CONTAINING_POINTS{1} =  POLYINFOloc.ELEMENTS_CONTAINING_xNEW   ;
%HISTORY.Bhred_interp{1} =  POLYINFOloc.BdomRED_interp   ;

DATA_ECM = DefaultField(DATA_ECM,'criterion_for_removing_weights','STANDARD' ); %  'ONLYweights'
DATA_ECM = DefaultField(DATA_ECM,'FORCE_POINTS_TO_REMAIN_INSIDE',1); %   Force points to remain inside of the domain
DATA_ECM = DefaultField(DATA_ECM,'THRESHOLD_NUMBER_OF_NEGATIVE_WEIGHTS',5); %  Negative points allowed at each iteration
DATA_ECM = DefaultField(DATA_ECM,'LINE_SEARCH_AVOID_NEGATIVE_WEIGTHS',0); %  cONTROL THE SIZE OF THE STEP ALONG THE DESCENT DIRECTION
DATA_ECM = DefaultField(DATA_ECM,'USE_LEAST_NORM_SOLUTION',0); %  Least-norm solution (pseudo-inverse)
DATA_ECM = DefaultField(DATA_ECM,'SECOND_STAGE_ITERATIONS',0); % Refine iterations, 2nd stage
DATA_ECM = DefaultField(DATA_ECM,'ONE_STAGE_ITERATIONS',0); % Refine iterations directly, in the first stage
DATA_ECM = DefaultField(DATA_ECM,'EXCLUDE_ELEMENT_TRANSITION_LIST',[]); % Freeze the positions of some points when they try to move (checked at each iteration)
DATA_ECM = DefaultField(DATA_ECM,'maxITER_allowed_residual_withoutDECREASE',10); %  Number of iterations allowed without the residual experiencing any decrease
DATA_ECM = DefaultField(DATA_ECM,'ELEMENTS_NOT_ALLOWED_TO_MOVE',[]); %  Similar to EXCLUDE_ELEMENT_TRANSITION_LIST,
% but imposed invariably at each iteration (so that the points cannot move even insided the element)
DATA_ECM = DefaultField(DATA_ECM,'EXCLUDE_ELEMENTS_TRIGGERED_NONCONVERGENCE_ITERATION',0); % When greater than 1, it automatically construct the
% list of elements EXCLUDE_ELEMENT_TRANSITION_LIST

if ~isempty(DATA_ECM.ELEMENTS_NOT_ALLOWED_TO_MOVE)
    ListElements = load(DATA_ECM.ELEMENTS_NOT_ALLOWED_TO_MOVE) ;
    DATA_ECM.ELEMENTS_NOT_ALLOWED_TO_MOVE = ListElements(:,1) ;
    PointsFixed = small2large(DATA_ECM.ELEMENTS_NOT_ALLOWED_TO_MOVE,VAR_SMOOTH_FE.ngausE) ;
    [~,~,DATA_ECM.FixedPointsAllIterations ]= intersect(PointsFixed,VAR_SMOOTH_FE.setPoints);
else
    DATA_ECM.FixedPointsAllIterations = [] ;
end



% to cross to other elements in the list

% List of design variables (DOFs )
% --------------------------------
[npoints, ~ ]= size(xOLD) ;
VARC.POINTSRpFIXED =  DATA_ECM.FixedPointsAllIterations; % These are the indexes of the points for which   position   are constrained, but weights aren't. These
% conditions doest not vary during the iterations
VARC.POINTSRp =VARC.POINTSRpFIXED ; % These are the indexes of the points for which   position   are constrained, but weights aren't. These conditions might
% vary during the iterations. At least, they must contain VARC.POINTSRpFIXED
VARC.POINTSl = 1:npoints; % These are the indexes of the points for which both weights and positions are unknowns
VARC.POINTSl(VARC.POINTSRp)  = [] ;


VARC.POINTSRpw = []; % These are the indexes of the points for which   weights and positions are constrained


% if DATA_ECM.ONE_STAGE_ITERATIONS>1
%     error('This option did not bring any advantage: set it == 0')
%     DATA_ECM.SECOND_STAGE_ITERATIONS  =  DATA_ECM.ONE_STAGE_ITERATIONS ;
%     [xSECONDS,wSECONDS,DATA_ECM,POLYINFO,VARC,HISTORY,ELEMENTS_TO_PLOT] = ControlPoints2ndSTAGE(xOLD,wOLD,b,DATA_ECM,VAR_SMOOTH_FE,VARC,HISTORY,ELEMENTS_TO_PLOT) ;
%
%
% else
[xFIRSTS,wFIRSTS,DATA_ECM,POLYINFO,VARC,HISTORY,ELEMENTS_TO_PLOT] = ControlPointsStandard(xOLD,wOLD,b,DATA_ECM,VAR_SMOOTH_FE,VARC,HISTORY,ELEMENTS_TO_PLOT) ;

if DATA_ECM.SECOND_STAGE_ITERATIONS >1 && length(find(wFIRSTS>0)) >1
    [xSECONDS,wSECONDS,DATA_ECM,POLYINFO,VARC,HISTORY,ELEMENTS_TO_PLOT] = ControlPoints2ndSTAGE(xFIRSTS,wFIRSTS,b,DATA_ECM,VAR_SMOOTH_FE,VARC,HISTORY,ELEMENTS_TO_PLOT) ;
end
% end


% Iterations with positive weights
% ----------------------------------
IterPositive= find(HISTORY.ISALLPOSITIVE ==1) ;
LastIteration = IterPositive(end) ;
% WE take this iteration as the final set of points
xNEW = HISTORY.POINTS{LastIteration} ;
wNEW = HISTORY.WEIGHTS{LastIteration} ;

HISTORY.POINTS =  HISTORY.POINTS(1:LastIteration)  ;
HISTORY.WEIGHTS =  HISTORY.WEIGHTS(1:LastIteration)  ;
HISTORY = DefaultField(HISTORY,'INDEXES_FIRST_STAGE',[]) ; 
HISTORY = DefaultField(HISTORY,'INDEXEX_SECOND_STAGE',[]) ; 

if  isempty(HISTORY.INDEXES_FIRST_STAGE)
    DATA_ECM.Include2ndStageIterations_PlotEvolutionWeights = 0  ; 
end

if DATA_ECM.Include2ndStageIterations_PlotEvolutionWeights ==1 
    IndSecondStageTake = LastIteration-length(HISTORY.INDEXES_FIRST_STAGE) ;
    if IndSecondStageTake <=0
        HISTORY.INDEXES_FIRST_STAGE = HISTORY.INDEXES_FIRST_STAGE(1:LastIteration) ;
        HISTORY.INDEXEX_SECOND_STAGE = {} ;
        HISTORY.POINTS_all =  HISTORY.POINTS_all(1:LastIteration)  ;
    HISTORY.WEIGHTS_all =  HISTORY.WEIGHTS_all(1:LastIteration)  ;
    else
        HISTORY.INDEXEX_SECOND_STAGE = HISTORY.INDEXEX_SECOND_STAGE(1:IndSecondStageTake) ;
        LastIndex = HISTORY.INDEXEX_SECOND_STAGE{end}(end) ; 
        HISTORY.POINTS_all =  HISTORY.POINTS_all(1:LastIndex)  ;
    HISTORY.WEIGHTS_all =  HISTORY.WEIGHTS_all(1:LastIndex)  ;
        
    end
else
     
    HISTORY.POINTS_all =  HISTORY.POINTS_all(1:LastIteration)  ;
    HISTORY.WEIGHTS_all =  HISTORY.WEIGHTS_all(1:LastIteration)  ;
end





HISTORY.ELEMENTS_CONTAINING_POINTS =  HISTORY.ELEMENTS_CONTAINING_POINTS(1:LastIteration)  ;
%HISTORY.Bhred_interp = HISTORY.Bhred_interp(1:LastIteration) ;

disp('---------------------------------')
disp('Summary')
disp('-------------------------------')
disp(['Reduction in first stage: from ',num2str(length(wFIRSTS)),' to ',num2str(length(find(wFIRSTS ~=0)))])

disp('--------------------------------------------------')
disp(['Final integration rule with m =',num2str(length(wNEW)),' POINTS  (of ',num2str(size(coorECM,1)),'). ',...
    '. Rank Basis = ',num2str(length(b))]);
disp('--------------------------------------------------') ;


disp('Integration error') ;

%disp(DATA_ECM.MSGPRINT{end})

%[PHIk_y,~, POLYINFO]= EvaluateBasisFunctionAtX_approx(xOLD, DATA_ECM.APPROX_FUN__DERI,VAR_SMOOTH_FE,POLYINFO)  ;

[PHIk_y,~,POLYINFO]=     EvaluateBasisFunctionALL(xNEW,DATA_ECM,VAR_SMOOTH_FE,POLYINFO) ;


bNEW = PHIk_y'*wNEW ;
errorINT = norm(bNEW-b)/norm(b)*100;
disp(['Error (%) =',num2str(errorINT)]) ;

%---------------
xINITIAL =coorECM; %(z,:) ;
wINITIAL = w ;
wGAUSS = wNEW ;
xGAUSS = xNEW ;
%ELEMENTS_xGAUSS= POLYINFO.ELEMENTS_CONTAINING_xNEW ;



 



% Nearest points
% ----------------
if size(coorECM,2) >1
    dt = DelaunayTri(COORg);
    NP = nearestNeighbor(dt, xNEW);
    
    xNEAR = COORg(NP,:); % Nearest point
else
    % 1D
    xNEAR = zeros(size(xNEW)) ;
    for iNEW = 1:length(xNEW)
        distALL = abs(COORg-xNEW(iNEW)) ;
        [~,indMIN] = min(distALL) ;
        xNEAR(iNEW) = COORg(indMIN) ;
    end
end

%[PHIk_y,~,~ ]= EvaluateBasisFunctionAtX_approx(xNEAR, DATA_ECM.APPROX_FUN__DERI,VAR_SMOOTH_FE,POLYINFO)  ;
[PHIk_y,~,~]=     EvaluateBasisFunctionALL(xNEAR,DATA_ECM,VAR_SMOOTH_FE,POLYINFO) ;


% Weights
wNEAR = PHIk_y'\b ;
bNEAR = PHIk_y'*wNEAR ;
errorINTnear = norm(bNEAR-b)/norm(b)*100;

 
disp(['Error nearest points (%) =',num2str(errorINTnear)]) ;
%disp(DATA_ECM.MSGPRINT{end}) ;
%%%

AUXVAR.PATHiniPOINTS = PATHiniPOINTS; 
AUXVAR.ELEMENTS_TO_PLOT = ELEMENTS_TO_PLOT; 
AUXVAR.HISTORY = HISTORY; 


 

