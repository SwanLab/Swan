function   [FACES_CN,rotDOM,angDOMall ]= ClassificationNodesMidSideCylin(CNmidside,COORmidside,CIRCLE,DATAIN)
if nargin ==0
    load('tmp3.mat')
end


% Converting cartesian coordinates to cylindrical coordinates
% ------------------------------------------------------------
xrel = COORmidside(:,1)-CIRCLE.CENTER(1) ;
zrel = COORmidside(:,3)-CIRCLE.CENTER(2) ;
yrel = COORmidside(:,2) ;
%%%
% Angle corresponding to each node
THETA = atan2(zrel,xrel) ;

%disp('This function is amenable to vectorization ... ')

FACES_CN = zeros(size(CNmidside)) ;
rotDOM = cell(1,size(CNmidside,1)) ; 
angDOMall = zeros(1,size(CNmidside,1)) ; 
for ielem = 1:size(CNmidside,1)
    % Checking rotation angle
    CNloc = CNmidside(ielem,:) ;
    ANGLE = THETA(CNloc,:) ;
    yLOC = yrel(CNloc,:) ;
    [dummy iface2] = min(yLOC) ;  % FACE 2 --> minimun y
    [dummy iface4] = max(yLOC) ; % FACE 4 --> maximum y
    % Now we have to choose faces 1 and 3
    INDICES = 1:4;
    INDICES([iface2 iface4]) = [] ;
    ANGLE_LOC = ANGLE(INDICES) ;
    
    [minANGLE indMIN_angle ]= min(ANGLE_LOC) ;
    [maxANGLE indMAX_angle ]= max(ANGLE_LOC) ;
    
    % By default, face 1 is the face with highest angle
    iface1 =  INDICES(indMAX_angle) ;
    iface3 =  INDICES(indMIN_angle) ;
    if maxANGLE-minANGLE >pi
        % In this case, the domain contains the region around 180 degrees
        % Hence
        iface3 =  INDICES(indMAX_angle) ;
        iface1 =  INDICES(indMIN_angle) ;
    end
    FACES_CN(ielem,:) = [iface1 iface2 iface3 iface4];
    
    % We can compute the rotation matrix here 
    % ----------------------------------------
    % Polar angle of face 1 
    AngleF1 = ANGLE(iface1) ; 
    % This angle is the angle subtended by the local e3 axis and the global
    % x-axis. The angle required for the matrix is obtained by subctracting
    % -pi/2 
  %  AngleDOMloc = AngleF1+pi/2; 
    % R3 = [cos(AngleDOMloc), 0 , sin(AngleDOMloc)
%      0            1        0 
%      -sin(AngleDOMloc)  0     cos(AngleDOMloc)] ; 

 
    
    rotDOM{ielem} = [sin(AngleF1), 0 , cos(AngleF1)
     0            1        0 
      -cos(AngleF1)  0     sin(AngleF1)] ;
    angDOMall(ielem) = AngleF1 ; 
 
    
    
    
    
end

end


%
%
% % Classification of midside nodes according to the faces they belong to
% % ----------------------------------------------------------------------
% xMIN = 1e20*ones(size(CNmidside,1),2) ;
% xMAX =-1e20*ones(size(CNmidside,1),2) ;
% for inodeE = 1:size(CNmidside,2)
%     COLUMN_loc = CNmidside(:,inodeE) ;
%     COOR_loc = COORmidside(COLUMN_loc,:) ;
%     for idim=1:2
%         xMIN(:,idim) = min([xMIN(:,idim), COOR_loc(:,idim)],[],2) ;
%         xMAX(:,idim) = max([xMAX(:,idim), COOR_loc(:,idim)],[],2) ;
%     end
% end
%
% FACES_CN = zeros(size(CNmidside)) ;
% for inodeE = 1:size(CNmidside,2)
%     COLUMN_loc = CNmidside(:,inodeE) ;
%     COOR_loc = COORmidside(COLUMN_loc,:) ;
%     IND =  (COOR_loc(:,1) == xMIN(:,1)) ;
%     FACES_CN(IND,inodeE) =  1 ;
%     IND =  (COOR_loc(:,1) == xMAX(:,1)) ;
%     FACES_CN(IND,inodeE) =  3 ;
%     IND =  (COOR_loc(:,2) == xMIN(:,2)) ;
%     FACES_CN(IND,inodeE) =  2 ;
%     IND =  (COOR_loc(:,2) == xMAX(:,2)) ;
%     FACES_CN(IND,inodeE) =  4 ;
% end