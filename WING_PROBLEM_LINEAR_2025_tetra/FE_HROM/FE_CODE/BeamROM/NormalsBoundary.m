function NORMALSv = NormalsBoundary(COOR,CNb)
% Compute normals of list CNb (connectivities)

if nargin == 0
    load('tmp.mat')
end

if isempty(CNb)
    error('You must generathe the Boundary Mesh !!!')
end

ndim = size(COOR,2) ;
NORMALSv = zeros(size(COOR,2),size(CNb,1));
Xloc1 = COOR(CNb(:,1),:)';
Xloc2 = COOR(CNb(:,2),:)';
t1 = Xloc2-Xloc1 ; % Tangential vector

if ndim == 3
    
    Xloc3 = COOR(CNb(:,3),:)';
    t2 = Xloc3-Xloc1 ; % Tangential vector
    n = cross(t1,t2) ;
    
else    
    n = zeros(size(t1)) ;
    n(1,:) = -t1(2,:) ;
    n(2,:) = t1(1,:) ;   
    
end

normN = sqrt(sum(n.^2,1)) ;
for idim = 1:size(n,1)
    NORMALSv(idim,:)  = n(idim,:)./normN ;
end

PRINT_LOCAL = 0; 

if PRINT_LOCAL ==1 
    figure(590)
    hold on 
    axis equal 
    plot(COOR(:,1),COOR(:,2),'r*') ;
    tiN = norm(t1(:,1)) ; 
    FACTOR = 3;
    NORMALSv = FACTOR*NORMALSv*tiN ; 
    
    for ielem = 1:size(CNb,1)
        DIFFCOOR = 0.5*(COOR(CNb(ielem,2),:) + COOR(CNb(ielem,1),:)) ; 
        
        xINI = DIFFCOOR' ; 
        xFIN = DIFFCOOR' + NORMALSv(:,ielem) ; 
        
        plot([xINI(1) xFIN(1)],[xINI(2), xFIN(2)],'r') ; 
    end
    
end




 end