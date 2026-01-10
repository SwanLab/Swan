function [xGAUSS,wGAUSS,DATALOC,ELEMENTS_xGAUSS,HISTORY,VAR_SMOOTH_FE,POLYINFO] ...
    = GenGauss2D3DfebasedCONTROL(COORg,w,DATALOC,coorECM,VAR_SMOOTH_FE)
%dbstop('6')
if nargin ==0
    load('tmp.mat')
   %  DATALOC.FORCE_POINTS_TO_REMAIN_INSIDE = 1 ; 
   %  DATALOC.SECOND_STAGE_ITERATIONS = 10 ; 
   %  DATALOC.ONE_STAGE_ITERATIONS = 0 ; 
end

b=   DATALOC.ExactIntegral; % PHI'*W ;
VAR_SMOOTH_FE.ExactIntegral = b; 

% Starting Cubature rule
% ----------------------
wOLD =w ;
xOLD = coorECM ;
PATHiniPOINTS = cell(size(xOLD,1),size(xOLD,1)) ;
DATALOC.mINI = size(xOLD,1) ;
 
for ipoint = 1:size(xOLD,1)
    PATHiniPOINTS{ipoint,1 } = xOLD(ipoint,:) ;
end
%INDEX_xNEW_respect_XOLD =  1:size(xOLD,1) ;
% At each iteration of the ensuing loop, a point of xOLD is removed
% We shall asign to the set of initial points the labels 1,2....npoints
% At each iteration, we have to ascertain which are the points that remains
% in the set
% For plotting affecgted elements
[~,~,POLYINFOloc]=     EvaluateBasisFunctionAtX_LARGE(xOLD,[],VAR_SMOOTH_FE,[]) ;
ELEMENTS_TO_PLOT =POLYINFOloc.ELEMENTS_CONTAINING_xNEW  ;


%NEW variables --> 24-Apr-2020 --> 
HISTORY.POINTS{1} = xOLD ; 
HISTORY.WEIGHTS{1} = wOLD ; 
HISTORY.ISALLPOSITIVE = 1  ;
HISTORY.ELEMENTS_CONTAINING_POINTS{1} =  POLYINFOloc.ELEMENTS_CONTAINING_xNEW   ;
%HISTORY.Bhred_interp{1} =  POLYINFOloc.BdomRED_interp   ;

DATALOC = DefaultField(DATALOC,'criterion_for_removing_weights','STANDARD' ); %  'ONLYweights'
DATALOC = DefaultField(DATALOC,'FORCE_POINTS_TO_REMAIN_INSIDE',1); %   Force points to remain inside of the domain 
DATALOC = DefaultField(DATALOC,'THRESHOLD_NUMBER_OF_NEGATIVE_WEIGHTS',5); %  Negative points allowed at each iteration
DATALOC = DefaultField(DATALOC,'LINE_SEARCH_AVOID_NEGATIVE_WEIGTHS',0); %  cONTROL THE SIZE OF THE STEP ALONG THE DESCENT DIRECTION
DATALOC = DefaultField(DATALOC,'USE_LEAST_NORM_SOLUTION',0); %  Least-norm solution (pseudo-inverse)
DATALOC = DefaultField(DATALOC,'SECOND_STAGE_ITERATIONS',0); % Refine iterations, 2nd stage
DATALOC = DefaultField(DATALOC,'ONE_STAGE_ITERATIONS',0); % Refine iterations directly, in the first stage
DATALOC = DefaultField(DATALOC,'EXCLUDE_ELEMENT_TRANSITION_LIST',[]); % Freeze the positions of some points when they try to move (checked at each iteration)
DATALOC = DefaultField(DATALOC,'maxITER_allowed_residual_withoutDECREASE',10); %  Number of iterations allowed without the residual experiencing any decrease
DATALOC = DefaultField(DATALOC,'ELEMENTS_NOT_ALLOWED_TO_MOVE',[]); %  Similar to EXCLUDE_ELEMENT_TRANSITION_LIST, 
% but imposed invariably at each iteration (so that the points cannot move even insided the element)
DATALOC = DefaultField(DATALOC,'EXCLUDE_ELEMENTS_TRIGGERED_NONCONVERGENCE_ITERATION',0); % When greater than 1, it automatically construct the
% list of elements EXCLUDE_ELEMENT_TRANSITION_LIST

if ~isempty(DATALOC.ELEMENTS_NOT_ALLOWED_TO_MOVE)
    ListElements = load(DATALOC.ELEMENTS_NOT_ALLOWED_TO_MOVE) ;
    DATALOC.ELEMENTS_NOT_ALLOWED_TO_MOVE = ListElements(:,1) ;
    PointsFixed = small2large(DATALOC.ELEMENTS_NOT_ALLOWED_TO_MOVE,VAR_SMOOTH_FE.ngausE) ; 
    [~,~,DATALOC.FixedPointsAllIterations ]= intersect(PointsFixed,VAR_SMOOTH_FE.setPoints);
else
    DATALOC.FixedPointsAllIterations = [] ; 
end

 

% to cross to other elements in the list

% List of design variables (DOFs )
% --------------------------------
[npoints, ~ ]= size(xOLD) ; 
VARC.POINTSRpFIXED =  DATALOC.FixedPointsAllIterations; % These are the indexes of the points for which   position   are constrained, but weights aren't. These
% conditions doest not vary during the iterations
VARC.POINTSRp =VARC.POINTSRpFIXED ; % These are the indexes of the points for which   position   are constrained, but weights aren't. These conditions might
% vary during the iterations. At least, they must contain VARC.POINTSRpFIXED 
VARC.POINTSl = 1:npoints; % These are the indexes of the points for which both weights and positions are unknowns 
VARC.POINTSl(VARC.POINTSRp)  = [] ; 


VARC.POINTSRpw = []; % These are the indexes of the points for which   weights and positions are constrained
 

if DATALOC.ONE_STAGE_ITERATIONS>1
    error('This option did not bring any advantage: set it == 0')
    DATALOC.SECOND_STAGE_ITERATIONS  =  DATALOC.ONE_STAGE_ITERATIONS ;
    [xSECONDS,wSECONDS,DATALOC,POLYINFO,VARC,HISTORY,ELEMENTS_TO_PLOT] = ControlPoints2ndSTAGE(xOLD,wOLD,b,DATALOC,VAR_SMOOTH_FE,VARC,HISTORY,ELEMENTS_TO_PLOT) ;


else
    [xFIRSTS,wFIRSTS,DATALOC,POLYINFO,VARC,HISTORY,ELEMENTS_TO_PLOT] = ControlPointsStandard(xOLD,wOLD,b,DATALOC,VAR_SMOOTH_FE,VARC,HISTORY,ELEMENTS_TO_PLOT) ;
    
    if DATALOC.SECOND_STAGE_ITERATIONS >1
        [xSECONDS,wSECONDS,DATALOC,POLYINFO,VARC,HISTORY,ELEMENTS_TO_PLOT] = ControlPoints2ndSTAGE(xFIRSTS,wFIRSTS,b,DATALOC,VAR_SMOOTH_FE,VARC,HISTORY,ELEMENTS_TO_PLOT) ;
    end
end


% Iterations with positive weights
% ----------------------------------
IterPositive= find(HISTORY.ISALLPOSITIVE ==1) ; 
LastIteration = IterPositive(end) ; 
% WE take this iteration as the final set of points 
xNEW = HISTORY.POINTS{LastIteration} ; 
wNEW = HISTORY.WEIGHTS{LastIteration} ; 

 HISTORY.POINTS =  HISTORY.POINTS(1:LastIteration)  ; 
 HISTORY.WEIGHTS =  HISTORY.WEIGHTS(1:LastIteration)  ; 
 HISTORY.ELEMENTS_CONTAINING_POINTS =  HISTORY.ELEMENTS_CONTAINING_POINTS(1:LastIteration)  ; 
%HISTORY.Bhred_interp = HISTORY.Bhred_interp(1:LastIteration) ; 



disp('--------------------------------------------------')
disp(['Final integration rule with m =',num2str(length(xNEW)),' POINTS  (of ',num2str(size(coorECM,1)),'). ',...
    '. Rank Basis = ',num2str(length(b))]);
disp('--------------------------------------------------') ;
 

disp('Integration error') ;

%disp(DATALOC.MSGPRINT{end})

%[PHIk_y,~, POLYINFO]= EvaluateBasisFunctionAtX_approx(xOLD, DATALOC.APPROX_FUN__DERI,VAR_SMOOTH_FE,POLYINFO)  ;

[PHIk_y,~,POLYINFO]=     EvaluateBasisFunctionAtX_LARGE(xNEW,[],VAR_SMOOTH_FE,POLYINFO) ; 

 
bNEW = PHIk_y'*wNEW ;
errorINT = norm(bNEW-b)/norm(b)*100;
disp(['Error (%) =',num2str(errorINT)]) ;

%---------------
xINITIAL =coorECM; %(z,:) ;
wINITIAL = w ;
wGAUSS = wNEW ;
xGAUSS = xNEW ;
ELEMENTS_xGAUSS= POLYINFO.ELEMENTS_CONTAINING_xNEW ; 


DATALOC.SHOW_FIGURES = 1;
if DATALOC.SHOW_FIGURES==1
    
    if size(coorECM,2) == 2
        PlotPathPointsGauss2D(ELEMENTS_TO_PLOT,VAR_SMOOTH_FE,PATHiniPOINTS,xINITIAL,xNEW) ; 
    else
        PlotPathPointsGauss3D(ELEMENTS_TO_PLOT,VAR_SMOOTH_FE,PATHiniPOINTS,xINITIAL,xNEW) ; 
    end
    
end

%dbstop('162')
% DATALOC = DefaultField(DATALOC,'HOLE',[]) ;
%
% if ~isempty(DATALOC.HOLE)
%     switch  DATALOC.HOLE.TYPE
%         case 'POLYGONAL'
%             %             POINTS_pol = DATA.DATAREMOVEPOINTS.HOLE.POINTS ;
%             %             POINTS_pol = [POINTS_pol;POINTS_pol(1,:)] ;
%             %             INPol = inpolygon(x,y,POINTS_pol(:,1),POINTS_pol(:,2)) ;
%             %             Xf(INPol,:) = 0 ;
%             POINTS_pol = DATALOC.HOLE.POINTS ;
%             POINTS_pol = [POINTS_pol;POINTS_pol(1,:)] ;
%             for ipoints = 1:size(POINTS_pol,1)-1
%                 plot(POINTS_pol(ipoints:ipoints+1,1),POINTS_pol(ipoints:ipoints+1,2),'b')
%             end
%
%         otherwise
%             error('Option not implemented')
%     end
% end





% Nearest points
% ----------------
dt = DelaunayTri(COORg);
NP = nearestNeighbor(dt, xNEW);

xNEAR = COORg(NP,:); % Nearest point
 
%[PHIk_y,~,~ ]= EvaluateBasisFunctionAtX_approx(xNEAR, DATALOC.APPROX_FUN__DERI,VAR_SMOOTH_FE,POLYINFO)  ;
[PHIk_y,~,~]=     EvaluateBasisFunctionAtX_LARGE(xNEAR,[],VAR_SMOOTH_FE,POLYINFO) ; 


 % Weights
wNEAR = PHIk_y'\b ;
bNEAR = PHIk_y'*wNEAR ;
errorINTnear = norm(bNEAR-b)/norm(b)*100;

%VAR_SMOOTH_FE.Bhred_interp = HISTORY.Bhred_interp{end} ;    % For checking error in interpolating the B-matrix
%VAR_SMOOTH_FE.Bhred_interp_INI = HISTORY.Bhred_interp{1} ;  

disp(['Error nearest points (%) =',num2str(errorINTnear)]) ;
%disp(DATALOC.MSGPRINT{end}) ;
%%%

if DATALOC.SHOW_FIGURES==1
    figure(16)
    hold on
    bar(sort(wNEW))
    ylabel('Weights (after optimization)')
    figure(17)
    hold on
    bar(sort(wINITIAL))
    ylabel('Weights (before optimization)')
    
end
%diary off
