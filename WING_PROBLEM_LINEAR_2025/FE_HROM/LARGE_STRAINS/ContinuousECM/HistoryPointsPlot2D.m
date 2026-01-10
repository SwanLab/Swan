function HistoryPointsPlot2D(AUXVAR,MESH,VAR_SMOOTH_FE)
% Patterned after
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/LARGE_STRAINS/ContinuousECM/HistoryPointsPlot1D.m

if nargin == 0
    load('tmp1.mat')
    close all
end

wALL = cell2mat(AUXVAR.HISTORY.WEIGHTS_all) ;
xALL =  AUXVAR.HISTORY.POINTS_all ;
%
COOR = MESH.COOR ;
xMIN = min(COOR(:,1)) ;
xMAX = max(COOR(:,1)) ;

yMIN = min(COOR(:,2)) ;
yMAX = max(COOR(:,2)) ;

zMAX  = max(max(wALL)) ;
AMP = AUXVAR.AmplificationFactorDomain;

zLIM = [0,zMAX+ AMP*(zMAX)] ;
%
% AMPL = 1.3 ;
% DATA.LIMITS = [AMPL*xMIN,AMPL*xMAX,AMPL*yMIN,AMPL*yMAX,zLIM] ;
DATA.LIMITS = [xMIN - AMP*(xMAX-xMIN),xMAX+ AMP*(xMAX-xMIN),yMIN- AMP*(yMAX-yMIN),yMAX+AMP*(yMAX-yMIN),zLIM] ;
%
AUXVAR = DefaultField(AUXVAR,'fhandle',[]) ; %  =f ;
if isempty(AUXVAR.fhandle)
    f=  figure ;   
else
    f = AUXVAR.fhandle ;
end

 
hold on
%ax = axes('Parent',f,'position',[0.13 0.39  0.77 0.54]);
hold on
iter = 1;
h = zeros(size(wALL,1),1) ;
hBASE = zeros(size(wALL,1),1) ;
AUXVAR = DefaultField(AUXVAR,'LineWidth',4) ;
ColorLOC = AUXVAR.ColorLines ; 

for iplot = 1:length(h)
    
    if AUXVAR.ShowLocationPointsBASE == 1; 
         hBASE(iplot) = plot3(xALL{iter}(iplot,1) ,xALL{iter}(iplot,2),0,...
        'Color',ColorLOC,'Marker','o','MarkerSize',6,'MarkerFaceColor',ColorLOC);
    else
        hBASE = [] ; 
    end
    
    
    h(iplot) = plot3(xALL{iter}(iplot,1)*ones(1,2) ,xALL{iter}(iplot,2)*ones(1,2),[0,wALL(iplot,iter)],...
        'Color',ColorLOC,'LineWidth',AUXVAR.LineWidth);
%     if iter == length(xALL)
%         plot3(xALL{iter}(iplot,1),xALL{iter}(iplot,2),0,'rx','MarkerSize',6);
%     end
end
DATA.h_cecm = h(1) ;
DATA.h_cecm_base = hBASE(1) ;

axis(DATA.LIMITS)  ;
axis manual
ASPEC_RATIO = AUXVAR.ASPECT_RATIO ;
daspect(ASPEC_RATIO) ;
disp(['Employed aspect ratio (plotting purposes):  daspect(',num2str(ASPEC_RATIO),')']);


xlabel('x')
ylabel('y')
zlabel('Weights')

hh = quadplot(MESH.CN,MESH.COOR(:,1),MESH.COOR(:,2)) ;
DATA.LEGEND_GAUSS = [] ;

AUXVAR = DefaultField(AUXVAR,'DATA_from_MAIN',[]) ;
AUXVAR.DATA_from_MAIN = DefaultField(AUXVAR.DATA_from_MAIN,'SHOW_GAUSS_POINTS_IN_ITERATIVE_PLOT',0) ;

AUXVAR.DATALOC = DefaultField(AUXVAR.DATALOC,'SHOW_GAUSS_POINTS_IN_ITERATIVE_PLOT',0) ; % = 0 ;
if AUXVAR.DATA_from_MAIN.SHOW_GAUSS_POINTS_IN_ITERATIVE_PLOT == 0
    AUXVAR.DATALOC.xGAUSS = [] ;
end


if ~isempty(AUXVAR.DATALOC.xGAUSS)
    xGAUSS=  AUXVAR.DATALOC.xGAUSS;
    wGAUSS=  AUXVAR.DATALOC.wGAUSS;
    hg = zeros(size(wGAUSS)) ;
    for iplot = 1:length(hg)
        hg(iplot) = plot3(xGAUSS(iplot,1)*ones(1,2) ,xGAUSS(iplot,2)*ones(1,2),[0,wGAUSS(iplot)],'Color',[1,0,0],'LineWidth',AUXVAR.LineWidth);
    end
    DATA.LEGEND_GAUSS = ['Gauss rule (',num2str(length(wGAUSS)),')'] ;
    DATA.h_gauss = hg(1) ;
end



%dd = get(h,'BarWidth') ;
% Plot in red the next point to be eliminated
if size(wALL,2) == 1
    
else 
iterNEXT = iter + 1;
IndNOW = find(wALL(:,iter) == 0) ;
npoints = size(wALL,1)-length(IndNOW) ;
IndNEXT = find(wALL(:,iterNEXT) == 0) ;
IndRemoveLOC = setdiff(IndNEXT,IndNOW) ;
xRED = xALL{iter}(IndRemoveLOC,:) ;
% 
%  AUXVAR.ColorToBeEliminated  = [0,0,0] ; 
%         AUXVAR.MarkerToBeEliminated  = 'o' ; 
%  AUXVAR.MarkerSizeToBeEliminated  = 6 ; 


hPOINT = plot3(xRED(1),xRED(2),0,AUXVAR.MarkerToBeEliminated,'MarkerSize',AUXVAR.MarkerSizeToBeEliminated,'Color',AUXVAR.ColorToBeEliminated,'LineWidth',3);

 
if AUXVAR.DATA_from_MAIN.HIDE_LEGEND == 0
    if ~isempty(DATA.LEGEND_GAUSS)
        DATA.legend = legend([hPOINT,DATA.h_cecm,DATA.h_gauss],{'Point to be eliminated','CECM rule',DATA.LEGEND_GAUSS});
    else
        DATA.legend = legend([hPOINT,DATA.h_cecm],{'Point to be eliminated','CECM rule'});
    end
    
else
    
    DATA.legend = [] ; 
    
end


htitle = title(['First stage (elimination): Number of points = ',num2str(npoints),' (of ',num2str(size(wALL,1)),')'])  ;
grid on

icluster = 1;

DATA.INDEXES_FIRST_STAGE  = AUXVAR.HISTORY.INDEXES_FIRST_STAGE ;
DATA.INDEXES_SECOND_STAGE  = AUXVAR.HISTORY.INDEXEX_SECOND_STAGE ;
AUXVAR.HISTORY = DefaultField(AUXVAR.HISTORY,'CONTROL_POINTS',[]) ;
DATA.CONTROL_POINTS = AUXVAR.HISTORY.CONTROL_POINTS ;
DATA.Include2ndStageIterations_PlotEvolutionWeights  = AUXVAR.DATALOC.Include2ndStageIterations_PlotEvolutionWeights ;



DATA.wALL = wALL  ;
DATA.xALL = xALL  ;

%npointsEND =  length(find(wALL(:,end) ~= 0)) ;
niter = size(wALL,2) ;% - npointsEND+1 ;
DATA.niter = niter;

DATA.h = h  ;
DATA.hBASE = hBASE; 
DATA.hPOINT = hPOINT;
DATA.htitle = htitle ;
guidata(gcf,DATA);


%
% SliderStep = [1/(niter-1),1/(niter-1)]  ;
%
% b = uicontrol('Parent',f,'Style','slider','Position',[81,20,419,23],...
%     'value',icluster, 'min',1, 'max',niter,'Callback','UpdateWeights2D','SliderStep',SliderStep,'Tag','slide_tag');

AUXVAR = DefaultField(AUXVAR,'DATA_from_MAIN',[]) ;
AUXVAR.DATA_from_MAIN = DefaultField(AUXVAR.DATA_from_MAIN,'MAKE_VIDEO_POINTS',0) ;


if   AUXVAR.DATA_from_MAIN.MAKE_VIDEO_POINTS ==1
    
    
    AUXVAR.DATA_from_MAIN = DefaultField(AUXVAR.DATA_from_MAIN,'NameVideo','animation.gif') ;
    DATA.NameVideo = AUXVAR.DATA_from_MAIN.NameVideo ;
    
    AUXVAR.DATA_from_MAIN = DefaultField(AUXVAR.DATA_from_MAIN,'DelayBetweenFrames',0.001) ;
    DATA.DelayBetweenFrames = AUXVAR.DATA_from_MAIN.DelayBetweenFrames ;
    
    AUXVAR = DefaultField(AUXVAR,'angle1',[]) ;
    AUXVAR = DefaultField(AUXVAR,'angle2',[]) ;
    
    DATA.anglesVIEW = [AUXVAR.angle1,AUXVAR.angle2] ;
    
    UpdateWeights2D_video(DATA,f) ;
else
    guidata(gcf,DATA);
    SliderStep = [1/(niter-1),1/(niter-1)]  ;
    b = uicontrol('Parent',f,'Style','slider','Position',[81,20,419,23],...
        'value',icluster, 'min',1, 'max',niter,'Callback','UpdateWeights2D','SliderStep',SliderStep,'Tag','slide_tag');
end

end

