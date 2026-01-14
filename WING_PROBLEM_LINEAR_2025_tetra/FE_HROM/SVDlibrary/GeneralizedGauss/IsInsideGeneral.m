function [inELEM,onELEM] = IsInsideGeneral(xLOC,COOR,CN,elemLOC,IND_POLYG)

ndim = size(COOR,2) ; 
CNloc = CN(elemLOC,:); % Nodes forming the element
if ndim == 2
    INDnodes =CNloc(IND_POLYG) ;  % Corner nodes (to define a polygon)
    COORelem = COOR(INDnodes,:) ; % Coordinates of the elemen
    [inELEM,onELEM] = inpolygon(xLOC(1),xLOC(2),COORelem(:,1),COORelem(:,2)) ;
else
    % Use triangulation
    % ------------------
    INDnodes =CNloc(IND_POLYG) ;  % Corner nodes
    INDnodes = INDnodes(1:end-1) ;
    COORelem = COOR(INDnodes,:) ; % Coordinates of the elemen
    TriLocal = delaunayTriangulation(COORelem);
    CHECKINSIDE = pointLocation(TriLocal,xLOC) ;
    if isnan(CHECKINSIDE)
        inELEM = 0 ; onELEM = 0 ;
    else
        inELEM = 1; onELEM = 0 ;
    end
    
    %           %
%           %
%           ic = incenter(TriLocal);
%           numtri = size(TriLocal,1);
%           trilabels = arrayfun(@(x) {sprintf('T%d', x)}, (1:numtri)');
%           Htl = text(ic(:,1), ic(:,2),ic(:,3), trilabels, 'FontWeight', ...
%               'bold', 'HorizontalAlignment', 'center', 'Color', ...
%               'blue');
%           plot3(xLOC(1),xLOC(2),xLOC(3),'x','MarkerSize',7)
%           hold off
%          
%          
end