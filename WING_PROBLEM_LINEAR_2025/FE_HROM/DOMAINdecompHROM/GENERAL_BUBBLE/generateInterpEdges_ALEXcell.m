function [COORquad, LINES,QUAD_DATA_REGION] = generateInterpEdges_ALEXcell(COOR,CN, p,QUADRANT_INDIVUAD)
 
% Initialize
if nargin == 0
    load('tmp1.mat')
end
COORquad = [];
LINES = cell(size(CN,1), 1);
nodeCount = 0;

for i = 1:length(LINES)
    
    nPoints = p + 1;

    idx1 = CN(i,1);  % Start point index
    idx2 = CN(i,2);      % End point index

    P1 = COOR(idx1, :);
    P2 = COOR(idx2, :);
    
    

    t = linspace(0, 1, nPoints)';
    edgePoints = (1 - t) * P1 + t * P2;
    
    if   i==1  || (i>1 && (i ~= size(CN,1)) && (CN(i-1,2) ~= CN(i,1))) 
        COORquad = [COORquad; edgePoints];
        LINES{i} = nodeCount + (1:nPoints);
        nodeCount = nodeCount + nPoints;
    elseif i>1 && (i ~= size(CN,1)) && CN(i-1,2) == CN(i,1) 
         COORquad = [COORquad; edgePoints(2:end,:)];
         LINES{i} = (nodeCount-1) + (1:nPoints);
         nodeCount = nodeCount + nPoints-1 ; 
    elseif  i == size(CN,1)
         COORquad = [COORquad; edgePoints(1:end-1,:)];
         
         LINES{i} = [nodeCount + (1:(nPoints-1)),1];
         
    end
    
    
   
end

% 
% 
% %%%% DETERMINE CENTROID QUADRANT REGIONS 
% QUAD_DATA_REGION= [] ; 
% for iregion = 1:size(QUADRANT_INDIVUAD)
%     
%    xmin = min(COOR(QUADRANT_INDIVUAD(iregion,:),1)) ; 
%    xmax = max(COOR(QUADRANT_INDIVUAD(iregion,:),1)) ; 
%    Lx = 0.5*(xmax-xmin) ; 
%    
%    ymin = min(COOR(QUADRANT_INDIVUAD(iregion,:),2)) ; 
%    ymax = max(COOR(QUADRANT_INDIVUAD(iregion,:),2)) ; 
%    Ly = 0.5*(ymax-ymin) ; 
%    
%    CENTROID = 0.5*[xmin+xmax,ymin+ymax] ; 
% %    
% %   %  QUAD_DATA_REGION(iregion) =[] ; 
% %      QUAD_DATA_REGION(iregion).CENTROID = CENTROID ; 
% %       QUAD_DATA_REGION(iregion).LIMS = [xmin ymin; xmax ymax] ;
% %        QUAD_DATA_REGION(iregion).LengthScale = [Lx,Ly] ; 
%      
% end
%  
