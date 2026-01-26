function COORquadNEW = COORheaFROMfile( CENTROID,DATAcommon) ;
% Nodes hexahedra element
% JAHO, 15-April-2024, Balmes 185, Barcelona 
%------------------------------
if nargin == 0
    load('tmp.mat')
end
 
if ~isempty(MESH_QUADRATIC_PARENT_DOMAIN)
    
    
    %% THESE ARE THE COORDINATES OF THE POINTS READ FROM THE TRAINING MESH FILE
%     COORgidPARENT = COORquad ;
%     ndim = size(COORgidPARENT,2);
%     
%     for idim = 1:ndim
%         COORgidPARENT(:,idim) = COORgidPARENT(:,idim) + CENTROID(idim) ;
%     end
    % LET US REFER THEM TO ITS PSEUDO-CENTROID (computed using the corners, not the midsides)
%     COORgidPARENT  = COORgidPARENT(1:end-1,:) ;
%     CENTparent = sum(COORgidPARENT(1:4,:),1)/4;
%     for idim = 1:ndim
%         COORgidPARENT(:,idim) = COORgidPARENT(:,idim) - CENTparent(idim) ;
%     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % NOW WE READ THE COORDINATES OF THE PARENT DOMAIN CREATED BY GID
%     DATAcommon = DefaultField(DATAcommon,'MESH_QUADRATIC_PARENT_DOMAIN',[]) ; 
%     if isempty(DATAcommon.MESH_QUADRATIC_PARENT_DOMAIN)
    MESHparent= ReadMeshFileStr(MESH_QUADRATIC_PARENT_DOMAIN,'READ_MATERIAL_COLUMN',0)  ;
%     else
%         MESHparent= ReadMeshFileStr(MESH_QUADRATIC_PARENT_DOMAIN,'READ_MATERIAL_COLUMN',1)  ;
%     end
     COORgid = MESHparent.COOR(MESHparent.CN,:) ;
%     if size(MESHparent.COOR,1) == 9
%         COORgid = COORgid(1:end-1,:) ;
%     end
    % LET US REFER THEM TO THE CENTROID
    
    
    
    
    CENT = CENTROID ; 
    ndim = size(CENTROID,2); 
    for idim = 1:ndim
        COORgid(:,idim) = COORgid(:,idim) - CENT(idim) ;
    end
    
    COORquadNEW = COORgid ;     
%     
%     k = dsearchn(COORgid,COORgidPARENT) ;
%     
%     DISTANCE= COORgidPARENT -COORgid(k,:)  ;
%     
%     distCORNERSnorm = norm(DISTANCE(1:4,:)) ;
    
%     if distCORNERSnorm > 1e-10
%         error('Incorrect parent domain')
%     else
%         % Therefore
%         COORquadNEW = COORgid(k,:)  ; ;
%         for idim = 1:ndim
%             COORquadNEW(:,idim) = COORquadNEW(:,idim) + CENTparent(idim) - CENTROID(idim) ;
%         end
%         COORquadNEW = [COORquadNEW; 0 0] ;
%         
%     end
    
else
    error('You must provide the input data file for the training boundary mesh')
end