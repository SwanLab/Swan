function HistoryPointsDECM3D(AUXVAR,MESH,VAR_SMOOTH_FE)
%--------------------------------------------------------------------------
% function HistoryPointsDECM3D(AUXVAR, MESH, VAR_SMOOTH_FE)
%
% PURPOSE:
%   Visualizes the evolution of weights and cubature points across 
%   iterations of the Continuous Empirical Cubature Method (CECM) in 3D 
%   domains. This includes a graphical representation of the active 
%   quadrature points and their relative weights at each iteration.
%
%   The weights are represented by the size of markers at each point 
%   location in 3D, offering a visual clue about their importance during 
%   the cubature optimization. Optionally, the evolution can be saved 
%   as a video using `MAKE_VIDEO_POINTS`.
%
% INPUTS:
%   - AUXVAR          : Structure containing plotting configurations and 
%                       history of the CECM process. Expected subfields:
%                         * DATAOUTdecm.HistoryPoints (cell array)
%                         * DATAOUTdecm.HistoryWeights (cell array)
%                         * HISTORY.WEIGHTS_all (full weights over all iterations)
%                         * DATA_from_MAIN.MAKE_VIDEO_POINTS (flag)
%                         * DATALOC (plot settings: marker sizes, etc.)
%
%   - MESH            : Mesh structure with at least the field:
%                         * COOR : Coordinates of the domain mesh nodes
%
%   - VAR_SMOOTH_FE   : Structure with smoothed data fields, including:
%                         * COORg : Global coordinates of Gauss or cubature points
%
% OUTPUT:
%   - None (plots directly to figure window). If `MAKE_VIDEO_POINTS` is
%     set to true, a `.gif` animation is generated showing the evolution.
%
% FUNCTIONAL OVERVIEW:
%   1. Loads DECM point histories and corresponding weights.
%   2. Reconstructs full history of selected points and aligns them with
%      global coordinates (`COORg`).
%   3. Computes a weight matrix `wALL`, assigning each point a marker
%      size according to its weight at a given iteration.
%   4. Displays domain mesh boundaries via `DomainContPLOT3D`.
%   5. Shows points with zero weight as transparent; weighted points as
%      markers sized proportionally.
%   6. If `MAKE_VIDEO_POINTS` is enabled, the point evolution is animated
%      using `UpdateWeights3Ddecm_video`, followed by final plot update.
%   7. Otherwise, a slider is created to manually navigate the iterations.
%
% DEPENDENCIES:
%   - DefaultField
%   - DomainContPLOT3D
%   - UpdateWeights3Ddecm_video
%   - HistoryPointsPlot3D
%   - VAR_SMOOTH_FE.COORg must correspond to all candidate points
%
% USAGE EXAMPLE:
%   HistoryPointsDECM3D(AUXVAR, MESH, VAR_SMOOTH_FE)
%
% NOTES:
%   - This function is intended for interactive or postprocessing analysis
%     of cubature rule construction and sparsification strategies.
%   - Axis scaling is preserved throughout iterations for clarity.
%
% SEE ALSO:
%   HistoryPointsALL3D, UpdateWeights3Ddecm_video, HistoryPointsPlot3D
%
% AUTHOR:
%   Joaquín A. Hernández, UPC-CIMNE, May 2025.
%--------------------------------------------------------------------------

% Patterned after
%/home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/LARGE_STRAINS/ContinuousECM/HistoryPointsPlot3D.m
if nargin == 0
    load('tmp1.mat')
    
   % AUXVAR.DATA_from_MAIN.MAKE_VIDEO_POINTS = 0 ;
    % AUXVAR.DATALOC.SHOW_GAUSS_POINTS_IN_ITERATIVE_PLOT = 0;
  % AUXVAR.DATALOC.SHOW_MESH_DOMAIN = 1; 
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


f=  figure ;
hold on

% Plot boundary of the domain (or the mesh)
 DomainContPLOT3D(AUXVAR,MESH,VAR_SMOOTH_FE) ; 



axis equal
hold on
iter = 1;
h = zeros(size(wALL,1),1) ;
% Magnitude of weights is represented as the size of the Marker
AUXVAR.DATALOC = DefaultField(AUXVAR.DATALOC,'PlotEvolutionWeights_MarkerSizeMin',5) ;
AUXVAR.DATALOC = DefaultField(AUXVAR.DATALOC,'PlotEvolutionWeights_MarkerSizeMax',100) ;
DATA.MarkerSizeMin = AUXVAR.DATALOC.PlotEvolutionWeights_MarkerSizeMin ;
DATA.MarkerSizeMax = AUXVAR.DATALOC.PlotEvolutionWeights_MarkerSizeMax ;


for iplot = 1:length(h)
    wPOINT = wALL(iplot,iter) ;
    
    if wPOINT == 0
        h(iplot) = plot3(xALL(iplot,1) ,xALL(iplot,2),xALL(iplot,3),'Color',[0,0,0],'Marker','none');
    else
        iplotREF = iplot;
        MarkerSizeLoc =  DATA.MarkerSizeMin +(DATA.MarkerSizeMax-DATA.MarkerSizeMin)/(DATA.wMAX-DATA.wMIN)*(wPOINT-DATA.wMIN) ;
        h(iplot) = plot3(xALL(iplot,1) ,xALL(iplot,2),xALL(iplot,3),'Color',[0,0,0],'Marker','.','MarkerSize',MarkerSizeLoc);
        
    end
end
DATA.h_cecm = h(iplotREF) ;

AUXVAR.CURRENT_AXIS = axis; 

xlabel('x')
ylabel('y')
zlabel('z')

axis(AUXVAR.CURRENT_AXIS );  

npoints = length(find(wALL(:,iter)) >0) ;

htitle = title(['DECM: error (%) =',num2str(DATA.errorDECM(iter)*100),';  Number of points = ',num2str(npoints),' (of ',num2str(size(wALL,1)),')'])  ;
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
    
    
    
    UpdateWeights3Ddecm_video(DATA,f) ;
    
    ok = menu('Press OK if you want to see the CECM reduction stages','OK','FINISH SIMULATION') ;
    
    
    [AUXVAR.angle1,AUXVAR.angle2 ] =view ;
    AUXVAR.fhandle = f ;
    
    if ok == 1
        
        for iplot = 1:length(DATA.h)
            delete(DATA.h(iplot))
        end
        delete(DATA.htitle) ;
        
        
        HistoryPointsPlot3D(AUXVAR,MESH,VAR_SMOOTH_FE) ;
    end
    
    
    
    
    
else
    guidata(gcf,DATA);
    SliderStep = [1/(niter-1),1/(niter-1)]  ;
    b = uicontrol('Parent',f,'Style','slider','Position',[81,20,419,23],...
        'value',icluster, 'min',1, 'max',niter,'Callback','UpdateWeights3D_DECM','SliderStep',SliderStep,'Tag','slide_tag');
end
