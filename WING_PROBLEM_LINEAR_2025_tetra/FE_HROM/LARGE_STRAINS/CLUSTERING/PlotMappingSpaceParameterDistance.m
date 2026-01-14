function DATAoffline = PlotMappingSpaceParameterDistance(DATAoffline,MATRIX_POINTS_SPACE_PARAMETER,INDEXES_TRAJECTORIES,SNAPdisp)

if nargin == 0
    load('tmp.mat')
end
DATAoffline = DefaultField(DATAoffline,'ifigureMANIFOLD',300) ; 
figure(DATAoffline.ifigureMANIFOLD)
hold on 
xlabel('Parameter 1')
ylabel('Parameter 2')
zlabel('Euclidean distance')

DISTANCES = [] ; 
PARAMETERS = [] ; 

DATAoffline = DefaultField(DATAoffline,'MarkerPLotMani','none')   ;  

if isempty(INDEXES_TRAJECTORIES)
    INDEXES_TRAJECTORIES = {1:size(SNAPdisp,2)};
end

for itraj = 1:length(INDEXES_TRAJECTORIES)
    
    INDEX_LOCAL = INDEXES_TRAJECTORIES{itraj} ; 
    SNAPloc = SNAPdisp(:,INDEX_LOCAL) ; 
    
    dSNAPloc = SNAPloc(:,2:end) - SNAPloc(:,1:end-1) ; 
    dSNAPloc = sqrt(sum(dSNAPloc.^2,1)) ; 
    
    PARAM = MATRIX_POINTS_SPACE_PARAMETER(INDEX_LOCAL(2:end),:) ;
    
    DISTANCES = [DISTANCES ; dSNAPloc'] ; 
    PARAMETERS = [PARAMETERS;PARAM ] ; 
    
    plot3(PARAM(:,1),PARAM(:,2),dSNAPloc,'LineWidth',2,'Marker',DATAoffline.MarkerPLotMani) ; 
    
end

% maxX = max(PARAMETERS,1) ; 
% minX = min(PARAMETERS,1) ;
% maxY = max(PARAMETERS,2) ; 
% minY = min(PARAMETERS,2) ;
% 
% ngrid = 50 ; 
% XX = minX:ngrid:maxX ; 
% YY = minY:ngrid:maxY ; 


% [xq,yq] = meshgrid(XX, YY);
% vq = griddata(PARAMETERS(:,1),PARAMETERS(:,2),DISTANCES,xq,yq);
% 
% mesh(xq,yq,vq)
% hold on
% plot3(x,y,v,'o')
 
