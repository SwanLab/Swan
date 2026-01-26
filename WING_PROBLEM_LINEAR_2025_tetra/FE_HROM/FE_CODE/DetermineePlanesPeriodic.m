function NODESpl =  DetermineePlanesPeriodic(COORbound,TOL,xmax,xmin,ymax,ymin,zmax,zmin,NODESbound)

if nargin == 0
    load('tmp.mat')
end
NODESpl = {} ; 
% 1) PLANES x = xmax
ref = xmax; idim = 1; 
%NODESpl.XMAX = find(abs(COORbound(:,idim)-ref)<TOL) ;
NODESpl{1,1} = find(abs(COORbound(:,idim)-ref)<TOL) ;; 
% 2) PLANE x = xmin
ref = xmin; idim = 1; 
NODESpl{1,2} = find(abs(COORbound(:,idim)-ref)<TOL) ;
% 3) PLANE y = ymax
ref = ymax; idim = 2; 
NODESpl{2,1} = find(abs(COORbound(:,idim)-ref)<TOL) ;
% 4) PLANE y = ymin
ref = ymin; idim = 2; 
NODESpl{2,2} = find(abs(COORbound(:,idim)-ref)<TOL) ;
% 5) PLANE z = zmax
ref = zmax; idim = 3; 
NODESpl{3,1} = find(abs(COORbound(:,idim)-ref)<TOL) ;
% 6) PLANE z = zmin
ref = zmin; idim = 3; 
NODESpl{3,2} = find(abs(COORbound(:,idim)-ref)<TOL) ;
%%% -----------------------------------------------
% Change from local numbering (boundary nodes) to global numbering

for i=1:3 
    for j=1:2
        NODESpl{i,j} = NODESbound(NODESpl{i,j}) ; 
    end
end 

% 
% % 1) PLANES x = xmax
% ref = xmax; idim = 1; 
% NODESpl.XMAX = find(abs(COORbound(:,idim)-ref)<TOL) ;
% % 2) PLANE x = xmin
% ref = xmin; idim = 1; 
% NODESpl.XMIN = find(abs(COORbound(:,idim)-ref)<TOL) ;
% % 3) PLANE y = ymax
% ref = ymax; idim = 2; 
% NODESpl.YMAX = find(abs(COORbound(:,idim)-ref)<TOL) ;
% % 4) PLANE y = ymin
% ref = ymin; idim = 2; 
% NODESpl.YMIN = find(abs(COORbound(:,idim)-ref)<TOL) ;
% % 5) PLANE z = zmax
% ref = zmax; idim = 3; 
% NODESpl.ZMAX = find(abs(COORbound(:,idim)-ref)<TOL) ;
% % 6) PLANE z = zmin
% ref = zmin; idim = 3; 
% NODESpl.ZMIN = find(abs(COORbound(:,idim)-ref)<TOL) ;
% %%% -----------------------------------------------