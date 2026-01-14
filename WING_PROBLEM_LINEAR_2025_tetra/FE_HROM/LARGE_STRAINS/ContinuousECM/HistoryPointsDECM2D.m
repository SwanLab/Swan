function HistoryPointsDECM2D(AUXVAR,MESH,VAR_SMOOTH_FE)
% Patterned after
%/home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/LARGE_STRAINS/ContinuousECM/HistoryPointsDECM3D.m
if nargin == 0
    load('tmp1.mat') 
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
xMIN = min(COOR(:,1)) ;
xMAX = max(COOR(:,1)) ;

yMIN = min(COOR(:,2)) ;
yMAX = max(COOR(:,2)) ;

zMAX  = max(max(wALL)) ;
AMPL  = 1.3; 
zLIM = [0,AMPL*zMAX] ;
%
 AMP = 0.2; 
DATA.LIMITS = [xMIN - AMP*(xMAX-xMIN),xMAX+ AMP*(xMAX-xMIN),yMIN- AMP*(yMAX-yMIN),yMAX+AMP*(yMAX-yMIN),zLIM] ;
iplotREF = 1; 

f=  figure ;
hold on

% % Plot boundary of the domain (or the mesh)
  
hh = quadplot(MESH.CN,MESH.COOR(:,1),MESH.COOR(:,2)) ;
DATA.LEGEND_GAUSS = [] ; 




axis equal
hold on
iter = 1;
h = zeros(size(wALL,1),1) ;
  
AUXVAR = DefaultField(AUXVAR,'LineWidth',4) ;

 
for iplot = 1:length(h)
    
    
    h(iplot) = plot3(xALL(iplot,1)*ones(1,2) ,xALL(iplot,2)*ones(1,2),[0,wALL(iplot,1)],...
        'Color',[0,0,0],'LineWidth',AUXVAR.LineWidth);
end

 
DATA.h_cecm = h(iplotREF) ;

axis(DATA.LIMITS)  ;
axis manual
ASPEC_RATIO = AUXVAR.ASPECT_RATIO ; 
 daspect(ASPEC_RATIO) ; 
 disp(['Employed aspect ratio (plotting purposes):  daspect(',num2str(ASPEC_RATIO),')']);
% 
% AUXVAR.CURRENT_AXIS = axis; 

xlabel('x')
ylabel('y')
zlabel('weight')

 

npoints = length(find(wALL(:,iter)) >0) ;

htitle = title(['DECM: error (%) =',num2str(DATA.errorDECM(iter)*100),';  Number of points = ',num2str(npoints),' (of ',num2str(length(HistoryPoints)),')'])  ;
grid on
icluster = 1;


DATA.wALL = wALL  ;
DATA.xALL = xALL  ;
niter = size(wALL,2) ;% - npointsEND+1 ;
DATA.niter = niter;
DATA.h = h  ;
%    DATA.hPOINT = hPOINT;
DATA.htitle = htitle ;


AUXVAR = DefaultField(AUXVAR,'DATA_from_MAIN',[]) ;
AUXVAR.DATA_from_MAIN = DefaultField(AUXVAR.DATA_from_MAIN,'MAKE_VIDEO_POINTS',0) ;


if   AUXVAR.DATA_from_MAIN.MAKE_VIDEO_POINTS ==1
    
    
    AUXVAR.DATA_from_MAIN = DefaultField(AUXVAR.DATA_from_MAIN,'NameVideo','animation.gif') ;
    DATA.NameVideo = AUXVAR.DATA_from_MAIN.NameVideo ;
    
    AUXVAR.DATA_from_MAIN = DefaultField(AUXVAR.DATA_from_MAIN,'DelayBetweenFrames',0.001) ;
    DATA.DelayBetweenFrames = AUXVAR.DATA_from_MAIN.DelayBetweenFrames ;
    
    DATA.NumberOfPointsTotal = length(HistoryPoints) ; 
    
    UpdateWeights2Ddecm_video(DATA,f) ;
    
    ok = menu('Press OK if you want to see the CECM reduction stages','OK','FINISH SIMULATION') ;
    
    
    [AUXVAR.angle1,AUXVAR.angle2 ] =view ;
    AUXVAR.fhandle = f ;
    
    if ok == 1
        
        for iplot = 1:length(DATA.h)
            delete(DATA.h(iplot))
        end
        delete(DATA.htitle) ;
        
        
        HistoryPointsPlot2D(AUXVAR,MESH,VAR_SMOOTH_FE) ;
    end
    
    
    
    
    
else
    guidata(gcf,DATA);
    SliderStep = [1/(niter-1),1/(niter-1)]  ;
    b = uicontrol('Parent',f,'Style','slider','Position',[81,20,419,23],...
        'value',icluster, 'min',1, 'max',niter,'Callback','UpdateWeights2D_DECM','SliderStep',SliderStep,'Tag','slide_tag');
end
