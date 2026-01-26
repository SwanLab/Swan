function HistoryPointsDECM1D(AUXVAR,MESH,VAR_SMOOTH_FE)
% Patterned after
if nargin == 0
    load('tmp1.mat')
    %  AUXVAR.DATA_from_MAIN.MAKE_VIDEO_POINTS = 0 ;
    % AUXVAR.DATALOC.SHOW_GAUSS_POINTS_IN_ITERATIVE_PLOT = 0;
    
end

% Build matrices of weights and points
HistoryPoints = AUXVAR.DATAOUTdecm.HistoryPoints ;
HistoryWeights = AUXVAR.DATAOUTdecm.HistoryWeights ;
DATA.errorDECM = AUXVAR.DATAOUTdecm.Error ;

%List of points
HistoryPoints_all = cellfun(@transpose,HistoryPoints,'UniformOutput',false) ;
HistoryPoints_all = unique(cell2mat(HistoryPoints_all(:)')) ;
zALL = HistoryPoints_all(:);
xALL = VAR_SMOOTH_FE.COORg(zALL,:) ; % Coordinates of all the DECM points
% Building wALL
niter = length(HistoryWeights) ;
wALL = zeros(length(zALL),niter) ;
for iter = 1:niter
    POINTloc = HistoryPoints{iter} ;
    wLOC = HistoryWeights{iter} ;
    [~,IB,IC] = intersect(zALL,POINTloc) ;
    wALL(IB,iter) = wLOC(IC) ;
end
%wMAX_1st  = max(max(wALL)) ;  % Maximum weight (scaling purpose)
wALL_2nd = cell2mat(AUXVAR.HISTORY.WEIGHTS_all) ; % Weights during iterations
DATA.wMAX = max(max(wALL_2nd)) ;  % Maximum weight (scaling purpose)
DATA.wMIN = 0 ;

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

iter = 1;
h = zeros(size(xALL,1)) ;
AUXVAR = DefaultField(AUXVAR,'LineWidth',4) ;
for iplot = 1:length(h)
    h(iplot) = plot([xALL(iplot),xALL(iplot)],[0,wALL(iplot)],'Color',[0,0,0],'LineWidth',AUXVAR.LineWidth);
end

axis(DATA.LIMITS)  ;
axis manual
npoints = length(find(wALL(:,iter)) >0) ;


htitle = title(['DECM: error (%) =',num2str(DATA.errorDECM(iter)*100),';  Number of points = ',num2str(npoints),' (of ',num2str(size(wALL,1)),')'])  ;


grid on

icluster = 1;



DATA.wALL = wALL  ;
DATA.xALL = xALL  ;

%npointsEND =  length(find(wALL(:,end) ~= 0)) ;
niter = size(wALL,2) ;% - npointsEND+1 ;
DATA.niter = niter;

DATA.h = h  ;
DATA.htitle = htitle ;
guidata(gcf,DATA);



if   AUXVAR.DATA_from_MAIN.MAKE_VIDEO_POINTS ==1
    
    
    AUXVAR.DATA_from_MAIN = DefaultField(AUXVAR.DATA_from_MAIN,'NameVideo','animation.gif') ;
    DATA.NameVideo = AUXVAR.DATA_from_MAIN.NameVideo ;
    
    AUXVAR.DATA_from_MAIN = DefaultField(AUXVAR.DATA_from_MAIN,'DelayBetweenFrames',0.01) ;
    DATA.DelayBetweenFrames = AUXVAR.DATA_from_MAIN.DelayBetweenFrames ;
    
    
    
    UpdateWeights1Ddecm_video(DATA,f) ;
    
    ok = menu('Press OK if you want to see the CECM reduction stages','OK','FINISH SIMULATION') ;
    
    
   % [AUXVAR.angle1,AUXVAR.angle2 ] =view ;
    AUXVAR.fhandle = f ;
    
    if ok == 1
        
        for iplot = 1:length(DATA.h)
            delete(DATA.h(iplot))
        end
        delete(DATA.htitle) ;
        
        AUXVAR.KeepOldGIF_file = 1; 
        HistoryPointsPlot1D(AUXVAR,MESH) ;
    end
    
    
    
    
    
else
    
    
    SliderStep = [1/(niter-1),1/(niter-1)]  ;
    
    b = uicontrol('Parent',f,'Style','slider','Position',[81,20,419,23],...
        'value',icluster, 'min',1, 'max',niter,'Callback','UpdateWeights1Ddecm','SliderStep',SliderStep,'Tag','slide_tag');
    
end

