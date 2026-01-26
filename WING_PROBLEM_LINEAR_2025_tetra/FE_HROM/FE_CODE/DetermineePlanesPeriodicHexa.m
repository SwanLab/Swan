function [NODESpl NORMALS]=  DetermineePlanesPeriodicHexa(COORbound,TOL,xmax,xmin,ymax,ymin,zmax,zmin,NODESbound)

% See HEXAG.jpg 
% NODESpl{i} --> NODES pertaining to planes i
% NORMALS{i} % Vector perpendicular to plane i so that COOR_i + NORMALS{i}  =COOR_i_SLAVE  
%

if nargin ==0
    load('tmp2.mat')
end

D = (xmax-xmin) ; 
NODESpl = cell(8,1) ;
NORMALS = NODESpl ;
% 1) PLANES x = xmax   (--> xM)
ref = xmax; idim = 1;
NODESpl{1} = find(abs(COORbound(:,idim)-ref)<TOL) ;
NORMALS{1} = D*[-1 0 0]' ;  % Vector perpendicular to plane 1 so that COOR_1 + NORMALS{1}  =COOR_1_SLAVE (PLANE 4)

% 2) PLANE x = xmin    (--> xm)
ref = xmin; idim = 1;
NODESpl{4} = find(abs(COORbound(:,idim)-ref)<TOL) ;
NORMALS{4} =  -NORMALS{1} ;
% 3) Plane xMAX-yMAX
% -------------------
PA = [((xmax+xmin)/2); ymax] ;  % Points pertaining to this plane
d = (xmax-xmin)/2*tand(30) ;
PB = [xmax  ; (PA(2) -d) ] ;
INPLANE = COORbound(:,2) - PA(2) - (PB(2) -PA(2))/(PB(1)-PA(1))*(COORbound(:,1) - PA(1)) ;
NODESpl{2} = find(abs(INPLANE)<TOL) ;
n = PA-PB ;
NORMALS{2} = D*[-n(2) n(1) 0]'/norm(n);  % Unit normal vector to plane 2

% 4) Plane xMIN-yMIN
% -------------------
PA = [((xmax+xmin)/2); ymin] ;  % Points pertaining to this plane
PB = [xmin  ; (PA(2) +d) ] ;
INPLANE = COORbound(:,2) - PA(2) - (PB(2) -PA(2))/(PB(1)-PA(1))*(COORbound(:,1) - PA(1)) ;
NODESpl{5} = find(abs(INPLANE)<TOL) ;
NORMALS{5} = - NORMALS{2} ;

% 5) Plane xMAX-yMIN
% -------------------
PA = [((xmax+xmin)/2); ymin] ;  % Points pertaining to this plane
PB = [xmax  ; (PA(2) +d) ] ;
INPLANE = COORbound(:,2) - PA(2) - (PB(2) -PA(2))/(PB(1)-PA(1))*(COORbound(:,1) - PA(1)) ;
NODESpl{6} = find(abs(INPLANE)<TOL) ;
n = PB-PA ;
NORMALS{6} = D*[-n(2) n(1) 0]'/norm(n);  % Unit normal vector to plane 2

% 6) Plane xMIN-yMAX
% -------------------
PA = [((xmax+xmin)/2); ymax] ;  % Points pertaining to this plane
PB = [xmin  ; (PA(2) -d) ] ;
INPLANE = COORbound(:,2) - PA(2) - (PB(2) -PA(2))/(PB(1)-PA(1))*(COORbound(:,1) - PA(1)) ;
NODESpl{3} = find(abs(INPLANE)<TOL) ;
NORMALS{3} = -NORMALS{6} ;

% 7) PLANES z = zmax
ref = zmax; idim = 3;
NODESpl{7} = find(abs(COORbound(:,idim)-ref)<TOL) ;
NORMALS{7} = [0 0 -1]'*(zmax-zmin) ;
% 8) PLANES z = zmin
ref = zmin; idim = 3;
NODESpl{8} = find(abs(COORbound(:,idim)-ref)<TOL) ;
NORMALS{8} =  -NORMALS{7} ;



for i=1:length(NODESpl)
    NODESpl{i} = NODESbound(NODESpl{i}) ;
end


end
