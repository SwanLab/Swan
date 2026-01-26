function HistoryPointsPlot1D(AUXVAR,MESH)
% Patterned after /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/CLUSTERING_HROM/10_HOMOG_3D_1T/WeightsPlotEvolution.m
% more specifically:
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/CLUSTERING_HROM/10_HOMOG_3D_1T/WeightPlot2D.m

if nargin == 0
    load('tmp.mat')
end

wALL = cell2mat(AUXVAR.HISTORY.WEIGHTS_all) ;
xALL = cell2mat(AUXVAR.HISTORY.POINTS_all) ;
%
COOR = MESH.COOR ;
xMIN = min(COOR) ;
xMAX = max(COOR) ;
xLIM = [xMIN,xMAX] ;
yMAX  = max(max(wALL)) ;
yLIM = [0,1.2*yMAX] ;
%
DATA.LIMITS = [xLIM,yLIM] ;




f=  figure ;
hold on


%ax = axes('Parent',f,'position',[0.13 0.39  0.77 0.54]);
hold on
iter = 1;
h = zeros(size(xALL,1)) ;
AUXVAR = DefaultField(AUXVAR,'LineWidth',4) ;
for iplot = 1:length(h)
    h(iplot) = plot([xALL(iplot,iter),xALL(iplot,iter)],[0,wALL(iplot,iter)],'Color',[0,0,0],'LineWidth',AUXVAR.LineWidth);
end

axis(DATA.LIMITS)  ;
axis manual

%dd = get(h,'BarWidth') ;
% Plot in red the next point to be eliminated
iterNEXT = iter + 1;
IndNOW = find(wALL(:,iter) == 0) ;
npoints = size(wALL,1)-length(IndNOW) ;
IndNEXT = find(wALL(:,iterNEXT) == 0) ;
IndRemoveLOC = setdiff(IndNEXT,IndNOW) ;
xRED = xALL(IndRemoveLOC,iter) ;
wRED = wALL(IndRemoveLOC,iter) ;
hPOINT = plot(xRED,0,'o','MarkerSize',10,'Color',[1,0,0]);

DATA.legend = legend(hPOINT,{'Point to be eliminated'})

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
DATA.hPOINT = hPOINT;
DATA.htitle = htitle ;
guidata(gcf,DATA);




if   AUXVAR.DATA_from_MAIN.MAKE_VIDEO_POINTS ==1
    
    
    AUXVAR.DATA_from_MAIN = DefaultField(AUXVAR.DATA_from_MAIN,'NameVideo','animation.gif') ;
    DATA.NameVideo = AUXVAR.DATA_from_MAIN.NameVideo ;
    
    AUXVAR.DATA_from_MAIN = DefaultField(AUXVAR.DATA_from_MAIN,'DelayBetweenFrames',0.001) ;
    DATA.DelayBetweenFrames = AUXVAR.DATA_from_MAIN.DelayBetweenFrames ;
    
    %     AUXVAR = DefaultField(AUXVAR,'angle1',[]) ;
    %     AUXVAR = DefaultField(AUXVAR,'angle2',[]) ;
    %
    %     DATA.anglesVIEW = [AUXVAR.angle1,AUXVAR.angle2] ;
    AUXVAR = DefaultField(AUXVAR,'KeepOldGIF_file',0) ; 
    DATA.KeepOldGIF_file = AUXVAR.KeepOldGIF_file ; 
    UpdateWeights1D_video(DATA,f) ;
else
    
    SliderStep = [1/(niter-1),1/(niter-1)]  ;
    
    b = uicontrol('Parent',f,'Style','slider','Position',[81,20,419,23],...
        'value',icluster, 'min',1, 'max',niter,'Callback','UpdateWeights1D','SliderStep',SliderStep,'Tag','slide_tag');
end