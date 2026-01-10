function [xGAUSS,wGAUSS,DATALOC,ELEMENTS_xGAUSS,HISTORY,VAR_SMOOTH_FE] ...
    = GenGauss2D3Dfebased(COORg,w,DATALOC,coorECM,VAR_SMOOTH_FE)
%dbstop('6')
if nargin ==0
    load('tmp1.mat')
end

DATALOC =  DefaultField(DATALOC,'MSGPRINT',{}) ;

b=   DATALOC.ExactIntegral; % PHI'*W ;
if ~isempty(VAR_SMOOTH_FE)
VAR_SMOOTH_FE.ExactIntegral = b; 
end

% Starting Cubature rule
% ----------------------
wOLD =w ;
xOLD = coorECM ;
PATHiniPOINTS = cell(size(xOLD,1),size(xOLD,1)) ;
DATALOC.mINI = size(xOLD,1) ;
CONVERGENCE = 1 ;
iter = 1;
DATALOC.iter = iter;
for ipoint = 1:size(xOLD,1)
    PATHiniPOINTS{ipoint,iter } = xOLD(ipoint,:) ;
end
INDEX_xNEW_respect_XOLD =  1:size(xOLD,1) ;
% At each iteration of the ensuing loop, a point of xOLD is removed
% We shall asign to the set of initial points the labels 1,2....npoints
% At each iteration, we have to ascertain which are the points that remains
% in the set
% For plotting affecgted elements
[~,~,POLYINFOloc]=     EvaluateBasisFunctionAtX_FEinterp(xOLD,[],VAR_SMOOTH_FE,[]) ;
ELEMENTS_TO_PLOT =POLYINFOloc.ELEMENTS_CONTAINING_xNEW  ;

%NEW variables --> 24-Apr-2020 --> 
HISTORY.POINTS{1} = xOLD ; 
HISTORY.WEIGHTS{1} = wOLD ; 
HISTORY.ISALLPOSITIVE = 1  ;
HISTORY.ELEMENTS_CONTAINING_POINTS{1} =  POLYINFOloc.ELEMENTS_CONTAINING_xNEW   ;
HISTORY.Bhred_interp{1} =  POLYINFOloc.BdomRED_interp   ;

DATALOC = DefaultField(DATALOC,'criterion_for_removing_weights','STANDARD' ); %  'ONLYweights'


% ------------------------------------------------------
while CONVERGENCE ==1 && length(xOLD)>1
    DATALOC.iter = iter ; 
    

    [xNEW,wNEW,CONVERGENCE,DATALOC,POLYINFO] = RemovingPointsAlg(xOLD,wOLD,b,DATALOC,VAR_SMOOTH_FE) ;
    if CONVERGENCE == 1
        
        ELEMENTS_TO_PLOT = [ELEMENTS_TO_PLOT;POLYINFO.ELEMENTS_CONTAINING_xNEW  ] ;
        HISTORY.POINTS{end+1} = xNEW ; 
        HISTORY.WEIGHTS{end+1} = wNEW ;
        HISTORY.ISALLPOSITIVE(end+1) = all(wNEW>0) ; 
        HISTORY.Bhred_interp{end+1} =  POLYINFO.BdomRED_interp   ;

        HISTORY.ELEMENTS_CONTAINING_POINTS{end+1} =  POLYINFO.ELEMENTS_CONTAINING_xNEW   ;
        iter = iter + 1 ;
        INDEX_xNEW_respect_XOLD(DATALOC.REMOVED_INDEX) = [] ;
        for  ipoint = 1:size(xNEW,1)
            PATHiniPOINTS{INDEX_xNEW_respect_XOLD(ipoint),iter } = xNEW(ipoint,:) ;
        end
        xOLD = xNEW ;
        wOLD = wNEW ;
    end
end

% Iterations with positive weights
% ----------------------------------
IterPositive= find(HISTORY.ISALLPOSITIVE ==1) ; 
LastIteration = IterPositive(end) ; 
% WE take this iteration as the final set of points 
xOLD = HISTORY.POINTS{LastIteration} ; 
wOLD = HISTORY.WEIGHTS{LastIteration} ; 

 HISTORY.POINTS =  HISTORY.POINTS(1:LastIteration)  ; 
 HISTORY.WEIGHTS =  HISTORY.WEIGHTS(1:LastIteration)  ; 
 HISTORY.ELEMENTS_CONTAINING_POINTS =  HISTORY.ELEMENTS_CONTAINING_POINTS(1:LastIteration)  ; 
HISTORY.Bhred_interp = HISTORY.Bhred_interp(1:LastIteration) ; 



DATALOC.MSGPRINT{end+1}='--------------------------------------------------' ;
disp(DATALOC.MSGPRINT{end})
DATALOC.MSGPRINT{end+1} = ['Final integration rule with m =',num2str(length(xOLD)),' POINTS  (of ',num2str(size(coorECM,1)),'). ',...
    '. Rank Basis = ',num2str(length(b))];
disp(DATALOC.MSGPRINT{end})
DATALOC.MSGPRINT{end+1}='--------------------------------------------------' ;
disp(DATALOC.MSGPRINT{end})


DATALOC.MSGPRINT{end+1}='Integration error' ;

disp(DATALOC.MSGPRINT{end})

[PHIk_y,~, POLYINFO]= EvaluateBasisFunctionAtX_approx(xOLD, DATALOC.APPROX_FUN__DERI,VAR_SMOOTH_FE,POLYINFO)  ;
 
bNEW = PHIk_y'*wOLD ;
errorINT = norm(bNEW-b)/norm(b)*100;
DATALOC.MSGPRINT{end+1}=['Error (%) =',num2str(errorINT)] ;
disp(DATALOC.MSGPRINT{end}) ;

%---------------
xINITIAL =coorECM; %(z,:) ;
wINITIAL = w ;
wGAUSS = wOLD ;
xGAUSS = xOLD ;
ELEMENTS_xGAUSS= POLYINFO.ELEMENTS_CONTAINING_xNEW ; 


DATALOC.SHOW_FIGURES = 1;
if DATALOC.SHOW_FIGURES==1
    
    if size(coorECM,2) == 2
        PlotPathPointsGauss2D(ELEMENTS_TO_PLOT,VAR_SMOOTH_FE,PATHiniPOINTS,xINITIAL,xOLD) ; 
    else
        PlotPathPointsGauss3D(ELEMENTS_TO_PLOT,VAR_SMOOTH_FE,PATHiniPOINTS,xINITIAL,xOLD) ; 
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
NP = nearestNeighbor(dt, xOLD);

xNEAR = COORg(NP,:); % Nearest point
 
[PHIk_y,~,~ ]= EvaluateBasisFunctionAtX_approx(xNEAR, DATALOC.APPROX_FUN__DERI,VAR_SMOOTH_FE,POLYINFO)  ;

 % Weights
wNEAR = PHIk_y'\b ;
bNEAR = PHIk_y'*wNEAR ;
errorINTnear = norm(bNEAR-b)/norm(b)*100;

VAR_SMOOTH_FE.Bhred_interp = HISTORY.Bhred_interp{end} ;    % For checking error in interpolating the B-matrix
VAR_SMOOTH_FE.Bhred_interp_INI = HISTORY.Bhred_interp{1} ;  

DATALOC.MSGPRINT{end+1}=['Error nearest points (%) =',num2str(errorINTnear)] ;
disp(DATALOC.MSGPRINT{end}) ;
%%%

if DATALOC.SHOW_FIGURES==1
    figure(16)
    hold on
    bar(sort(wOLD))
    ylabel('Weights (after optimization)')
    figure(17)
    hold on
    bar(sort(wINITIAL))
    ylabel('Weights (before optimization)')
    
end
%diary off
