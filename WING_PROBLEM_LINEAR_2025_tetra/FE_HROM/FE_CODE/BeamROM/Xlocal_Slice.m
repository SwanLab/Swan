function [xSIGN] = Xlocal_Slice(PROP,rr1,MESH1D,INDEXelem)

if nargin == 0
    load('tmp4.mat')
end
CHANGE_ORDER = 0 ; 
xSIGN = 1;
% PROP = DefaultField(PROP,'xLOCAL',1) ; %  Direction local x
% if ~isempty(PROP.xLOCAL)
%     xSIGN = PROP.xLOCAL ;
% else
%     xSIGN = 1;
% end

% Firstly, we check that the 3D slices has been defined so that face 2 is
% the x_f2 > x_f1
if rr1(1) < 0
    error('Face 2 should have higher x-coordinate')
end
xSIGN =1 ; 

% Now we check whether it is a curved element or not
if abs(norm(rr1(2:end))/rr1(1)) > 1e-6
    % Curved domain 
    % -----------------------
    CN = MESH1D.CN(INDEXelem,:) ;
    NODES = unique(CN(:)) ;  % Nodes of this portion of the 1D structure
    xNODES = MESH1D.COOR(NODES,:) ; % coordinates 
      [RADIUS,CENTER] = fit_circle_through_3_points(xNODES(1:3,1:2)) ; 
      
   % Now we take one element 
   CNloc = MESH1D.CN(INDEXelem(1),:) ;
   % Initial and final nodes, according to the current matrix CN
   nodoINI = CNloc(:,1) ; 
   nodoFIN = CNloc(:,2) ; 
   % Coordinates
   xNODES = MESH1D.COOR([nodoINI; nodoFIN],:) ;
   % Polar coordinates (angle)
   xrel = xNODES(:,1)-CENTER(1) ;
   yrel = xNODES(:,2)-CENTER(2) ;
   THETHA_ini  = atan2(yrel(1),xrel(1));
   THETHA_fin  = atan2(yrel(2),xrel(2));
     % Function atan2 returns angles ranging between -pi and +pi
   % The connectivities should be such that the sequence of nodes is made
   % clockwise, that is, the angle should diminish (except for the transition from -pi to pi)
   INC_angle = THETHA_fin-THETHA_ini ; 
    
   if abs(INC_angle) < 2*pi
       if INC_angle >0 
           CHANGE_ORDER = 1; 
       end
   else
       % Case in which THETA_1 = pi - incANG,  THETA_2 = -pi + incANG
       if INC_angle <0 
           CHANGE_ORDER = 1; 
       end 
   end
  
   
     
else
    % Straigh element. Check that
       CN = MESH1D.CN(INDEXelem,:) ;
    xFIN = MESH1D.COOR(CN(1,2),:) ;
    xINI = MESH1D.COOR(CN(1,1),:) ;
    
    rLOC = xFIN-xINI ;
    if rLOC(1) < 0
        % Change the order of the connectivity matrix
        CHANGE_ORDER =1;
        
    end
    
end

if CHANGE_ORDER == 1
   xSIGN = -1 ; 
end

