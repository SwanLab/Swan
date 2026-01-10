function [NODESpl NORMALS]=  DetermineePlanesPeriodicHexaBCS(COORbound,TOL,xmax,xmin,ymax,ymin,zmax,zmin,NODESbound)

% See 
% NODESpl{i} --> NODES pertaining to planes i
% NORMALS{i} % Vector perpendicular to plane i so that COOR_i + NORMALS{i}  =COOR_i_SLAVE
%  See HEXAG_2.jpg

if nargin ==0
    load('tmp2.mat')
end

if  ~isempty(zmax)
    
    D = (xmax-xmin) ;
    NODESpl = cell(6,1) ;
    NORMALS = NODESpl ;
    % 1) PLANES x = xmin   
    ref = xmin; idim = 1;
    NODESpl{1} = find(abs(COORbound(:,idim)-ref)<TOL) ;
    NORMALS{1} = D*[1 0 0]' ;  % Vector perpendicular to plane 1 so that COOR_1 + NORMALS{1}  =COOR_1_SLAVE (PLANE 4)
    
    % 2) PLANE x = xmax    
    ref = xmax; idim = 1;
    NODESpl{3} = find(abs(COORbound(:,idim)-ref)<TOL) ;
    NORMALS{3} =  -NORMALS{1} ;
    % 3) Plane ymin
    % -------------------
    ref = ymin; idim = 2;
    NODESpl{2} = find(abs(COORbound(:,idim)-ref)<TOL) ;
    NORMALS{2} = (ymax-ymin)*[0 1 0]' ;  % Vector perpendicular to plane 1 so that COOR_1 + NORMALS{1}  =COOR_1_SLAVE (PLANE 4)
    
    % 4) PLANE  ymax
    ref = ymax; idim = 2;
    NODESpl{4} = find(abs(COORbound(:,idim)-ref)<TOL) ;
    NORMALS{4} = -(ymax-ymin)*[0 1 0]' ;  % Vector perpendicular to plane 1 so that COOR_1 + NORMALS{1}  =COOR_1_SLAVE (PLANE 4)
    
    
    % ) PLANES z = zmax
    ref = zmax; idim = 3;
    NODESpl{5} = find(abs(COORbound(:,idim)-ref)<TOL) ;
    NORMALS{5} = [0 0 -1]'*(zmax-zmin) ;
    % ) PLANES z = zmin
    ref = zmin; idim = 3;
    NODESpl{6} = find(abs(COORbound(:,idim)-ref)<TOL) ;
    NORMALS{6} =  -NORMALS{5} ;
    
    
else
    
    D = (xmax-xmin) ;
    NODESpl = cell(4,1) ;
    NORMALS = NODESpl ;
    % 1) PLANES x = xmin   (--> xM)
    ref = xmin; idim = 1;
    NODESpl{1} = find(abs(COORbound(:,idim)-ref)<TOL) ;
    NORMALS{1} = D*[1 0 ]' ;  % Vector perpendicular to plane 1 so that COOR_1 + NORMALS{1}  =COOR_1_SLAVE (PLANE 4)
    
    % 2) PLANE x = xmin    (--> xm)
    ref = xmax; idim = 1;
    NODESpl{3} = find(abs(COORbound(:,idim)-ref)<TOL) ;
    NORMALS{3} =  -NORMALS{1} ;
    % 3) Plane ymin
    % -------------------
    ref = ymin; idim = 2;
    NODESpl{2} = find(abs(COORbound(:,idim)-ref)<TOL) ;
    NORMALS{2} = (ymax-ymin)*[0 1 ]' ;  % Vector perpendicular to plane 1 so that COOR_1 + NORMALS{1}  =COOR_1_SLAVE (PLANE 4)
    
    % 2) PLANE  ymax
    ref = ymax; idim = 2;
    NODESpl{4} = find(abs(COORbound(:,idim)-ref)<TOL) ;
    NORMALS{4} = -(ymax-ymin)*[0 1 ]' ;  % Vector perpendicular to plane 1 so that COOR_1 + NORMALS{1}  =COOR_1_SLAVE (PLANE 4)
    
    
    
    
    
end



for i=1:length(NODESpl)
    NODESpl{i} = NODESbound(NODESpl{i}) ;
end


end
