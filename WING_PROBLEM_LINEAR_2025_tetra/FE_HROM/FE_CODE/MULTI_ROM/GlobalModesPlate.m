function Q = GlobalModesPlate(COOR,BNODES,CENTROID)

if nargin ==2
    CENTROID = [0,0,0]' ; 
end

CENTROID = CENTROID(:) ;

ndim = size(COOR,2); 
COORb = COOR(BNODES,:) ; 
COORb = bsxfun(@minus,COORb',CENTROID)';


%  Global modes associated to the domain ( = 20) 
nnodeb = size(COORb,1) ; 
ndof = nnodeb*ndim ; 
nmodes = 20 ; 
Q = zeros(ndof,nmodes) ; 
% Modes 1 to 8  
% --------------
idim = 1; 
Mone = [ones(nnodeb,1),COORb(:,1),COORb(:,2),COORb(:,3),COORb(:,1).*COORb(:,2),COORb(:,3).*COORb(:,2),...
     COORb(:,1).*COORb(:,3),COORb(:,1).*COORb(:,2).*COORb(:,3)] ;  

Q(idim:ndim:end,1:8) = Mone ; 
% Modes 9 to 16  
% --------------
idim = 2; 
Q(idim:ndim:end,9:16) = Mone ; 
% Modes 17 to 20 
idim = 3; 
Mthree = [ones(nnodeb,1),COORb(:,1),COORb(:,2),COORb(:,1).*COORb(:,2)] ; 
Q(idim:ndim:end,17:20) = Mthree ; 