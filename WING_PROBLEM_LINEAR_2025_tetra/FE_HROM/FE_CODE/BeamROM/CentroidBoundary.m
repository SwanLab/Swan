function CENTROIDS = CentroidBoundary(COOR,CNb)
% Compute centroid of list CNb (connectivities)

% if isempty(CNb)
%     error('You must generathe the Boundary Mesh !!!')
% end

CENTROIDS = zeros(size(CNb,1),size(COOR,2));

nnodeE = size(CNb,2) ; 
ndim = size(COOR,2) ; 
for inode = 1:nnodeE 
    for idim = 1:ndim
    CENTROIDS(:,idim) = CENTROIDS(:,idim)  + COOR(CNb(:,inode),idim)  ; 
    end
end

CENTROIDS = CENTROIDS/nnodeE ; 
% 
% Xloc1 = COOR(CNb(:,1),:)';
% Xloc2 = COOR(CNb(:,2),:)';
% Xloc3 = COOR(CNb(:,3),:)';
% t1 = Xloc2-Xloc1 ; % Tangential vector
% t2 = Xloc3-Xloc1 ; % Tangential vector
% n = cross(t1,t2) ;
% normN = sqrt(sum(n.^2,1)) ;
% for idim = 1:size(n,1)
%     NORMALSv(idim,:)  = n(idim,:)./normN ;
% end