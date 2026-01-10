function   [FACES_CN  ] = ClassificationNodesVertCylin(CNvert,COORvert,CIRCLE) ;

% Classification of vertices nodes
% AC = 1, AD = 2, BC=3 , BC=4
if nargin == 0
    load('tmp.mat')
end


% Converting cartesian coordinates to cylindrical coordinates
% (of one single element)
% ------------------------------------------------------------
 xrel = COORvert(:,1)-CIRCLE.CENTER(1) ;
zrel = COORvert(:,3)-CIRCLE.CENTER(2) ;
yrel = COORvert(:,2) ;
%%%
% Angle corresponding to each node
THETA = atan2(zrel,xrel) ;


% Checking polar  angle
ielem = 1; 
CNloc = CNvert(ielem,:) ;
ANGLE = THETA(CNloc,:) ;
yLOC = yrel(CNloc,:) ;

% Either Nodes 1 or 4 --->  yMAX
[ymax ] = max(yLOC) ;  % FACE 2 --> minimun y
indYMAX= find(yLOC == ymax) ; 
% Node 1 --> Higher angle, Node 4 --> Lowe angle
ANGLE_LOC = ANGLE(indYMAX) ;
[minANGLE indMIN_angle ]= min(ANGLE_LOC) ;
[maxANGLE indMAX_angle ]= max(ANGLE_LOC) ;
% By default, face 1 is the face with highest angle
icorner1 =  indYMAX(indMAX_angle) ;
icorner4 =  indYMAX(indMIN_angle) ;
if maxANGLE-minANGLE >pi
    % In this case, the domain contains the region around 180 degrees
    % Hence
    icorner4 =  indYMAX(indMAX_angle) ;
    icorner1 =  indYMAX(indMIN_angle) ;
end


% Either Nodes 2 or 3 --->  yMAX
[ymin ] = min(yLOC) ;  % FACE 2 --> minimun y
indYMIN = find(yLOC == ymin) ; 
% Node 1 --> Higher angle, Node 2 --> Lower angle , node 3
ANGLE_LOC = ANGLE(indYMIN) ;
[minANGLE indMIN_angle ]= min(ANGLE_LOC) ;
[maxANGLE indMAX_angle ]= max(ANGLE_LOC) ;
% By default, face 1 is the face with highest angle
icorner2 =  indYMIN(indMAX_angle) ;
icorner3 =  indYMIN(indMIN_angle) ;
if maxANGLE-minANGLE >pi
    % In this case, the domain contains the region around 180 degrees
    % Hence
    icorner3 =  indYMIN(indMAX_angle) ;
    icorner2 =  indYMIN(indMIN_angle) ;
end
 FACES_CN =[icorner1,icorner2,icorner3,icorner4] ; 

% 
% 
% % Now we have to choose faces 1 and 3
% INDICES = 1:4;
% INDICES([iface2 iface4]) = [] ;
% ANGLE_LOC = ANGLE(INDICES) ;
% 
% [minANGLE indMIN_angle ]= min(ANGLE_LOC) ;
% [maxANGLE indMAX_angle ]= max(ANGLE_LOC) ;
% 
% % By default, face 1 is the face with highest angle
% iface1 =  INDICES(indMAX_angle) ;
% iface3 =  INDICES(indMIN_angle) ;
% if maxANGLE-minANGLE >pi
%     % In this case, the domain contains the region around 180 degrees
%     % Hence
%     iface3 =  INDICES(indMAX_angle) ;
%     iface1 =  INDICES(indMIN_angle) ;
% end
% FACES_CN(ielem,[iface1 iface2 iface3 iface4]) = [1,2,3,4];
% 
% 
% 
% %
% %
% % % Classification of midside nodes according to the faces they belong to
% % % ----------------------------------------------------------------------
% % xMIN = 1e20*ones(size(CNmidside,1),2) ;
% % xMAX =-1e20*ones(size(CNmidside,1),2) ;
% % for inodeE = 1:size(CNmidside,2)
% %     COLUMN_loc = CNmidside(:,inodeE) ;
% %     COOR_loc = COORmidside(COLUMN_loc,:) ;
% %     for idim=1:2
% %         xMIN(:,idim) = min([xMIN(:,idim), COOR_loc(:,idim)],[],2) ;
% %         xMAX(:,idim) = max([xMAX(:,idim), COOR_loc(:,idim)],[],2) ;
% %     end
% % end
% %
% % FACES_CN = zeros(size(CNmidside)) ;
% % for inodeE = 1:size(CNmidside,2)
% %     COLUMN_loc = CNmidside(:,inodeE) ;
% %     COOR_loc = COORmidside(COLUMN_loc,:) ;
% %     IND =  (COOR_loc(:,1) == xMIN(:,1)) ;
% %     FACES_CN(IND,inodeE) =  1 ;
% %     IND =  (COOR_loc(:,1) == xMAX(:,1)) ;
% %     FACES_CN(IND,inodeE) =  3 ;
% %     IND =  (COOR_loc(:,2) == xMIN(:,2)) ;
% %     FACES_CN(IND,inodeE) =  2 ;
% %     IND =  (COOR_loc(:,2) == xMAX(:,2)) ;
% %     FACES_CN(IND,inodeE) =  4 ;
% % end
% 
% %
% % % Converting cartesian coordinates to cylindrical coordinates
% % % of just one single element
% % % ------------------------------------------------------------
% % ielem
% % xrel = COORmidside(:,1)-CIRCLE.CENTER(1) ;
% % zrel = COORmidside(:,3)-CIRCLE.CENTER(2) ;
% % yrel = COORmidside(:,2) ;
% % %%%
% % % Angle corresponding to each node
% % THETA = atan2(zrel,xrel) ;
% %
% %
% % FACES_CN = zeros(size(CNvert)) ;
% % for inodeE = 1:size(CNvert,2)
% %     COLUMN_loc = CNvert(:,inodeE) ;
% %     COOR_loc = COORvert(COLUMN_loc,:) ;
% %     IND_A =  (COOR_loc(:,1) == xMIN(:,1)) ;
% %     IND_B =  (COOR_loc(:,1) == xMAX(:,1)) ;
% %     IND_C =  (COOR_loc(:,2) == xMIN(:,2)) ;
% %     IND_D =  (COOR_loc(:,2) == xMAX(:,2)) ;
% %    % FACES_CN(IND,inodeE) =  4 ;
% %    %%%%% NODE 4-1
% %    IND = IND_A & IND_D ;
% %     FACES_CN(IND,inodeE) = 1 ;
% %      %%%%% NODE AC
% %    IND = IND_A & IND_C ;
% %     FACES_CN(IND,inodeE) = 2 ;
% %        %%%%% NODE BC
% %    IND = IND_B & IND_C ;
% %     FACES_CN(IND,inodeE) = 3 ;
% %         %%%%% NODE BD
% %    IND = IND_B & IND_D ;
% %     FACES_CN(IND,inodeE) = 4 ;
% %
% % end