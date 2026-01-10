function COORquadNEW = MidsideCoordintesQQ(COORquad,MESH_QUADRATIC_PARENT_DOMAIN,CENTROID) ;
% Recalculation of midside node coordinates, quadratic quadrilateral
% element
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/03_HETEROG/03_CURVED.mlx
% JAHO, 18-March-2023
%------------------------------
if nargin == 0
    load('tmp.mat')
end
COORquadNEW = COORquad;
if ~isempty(MESH_QUADRATIC_PARENT_DOMAIN)
    
    
    %% THESE ARE THE COORDINATES OF THE POINTS READ FROM THE TRAINING MESH FILE
    COORgidPARENT = COORquad ;
    ndim = size(COORgidPARENT,2);
    
    for idim = 1:ndim
        COORgidPARENT(:,idim) = COORgidPARENT(:,idim) + CENTROID(idim) ;
    end
    % LET US REFER THEM TO ITS PSEUDO-CENTROID (computed using the corners, not the midsides)
    COORgidPARENT  = COORgidPARENT(1:end-1,:) ;
    CENTparent = sum(COORgidPARENT(1:4,:),1)/4;
    for idim = 1:ndim
        COORgidPARENT(:,idim) = COORgidPARENT(:,idim) - CENTparent(idim) ;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % NOW WE READ THE COORDINATES OF THE PARENT DOMAIN CREATED BY GID
    MESHparent= ReadMeshFileStr(MESH_QUADRATIC_PARENT_DOMAIN,'READ_MATERIAL_COLUMN',0)  ;
    
    COORgid = MESHparent.COOR(MESHparent.CN,:) ;
    if size(MESHparent.COOR,1) == 9
        COORgid = COORgid(1:end-1,:) ;
    end
    % LET US REFER THEM TO ITS PSEUDO-CENTROID
    CENT = sum(COORgid(1:4,:),1)/4;
    for idim = 1:ndim
        COORgid(:,idim) = COORgid(:,idim) - CENT(idim) ;
    end
    
    k = dsearchn(COORgid,COORgidPARENT) ;
    
    DISTANCE= COORgidPARENT -COORgid(k,:)  ;
    
    distCORNERSnorm = norm(DISTANCE(1:4,:)) ;
    
    if distCORNERSnorm > 1e-8
        error('Incorrect parent domain')
    else
        % Therefore
        COORquadNEW = COORgid(k,:)  ; ;
        for idim = 1:ndim
            COORquadNEW(:,idim) = COORquadNEW(:,idim) + CENTparent(idim) - CENTROID(idim) ;
        end
        COORquadNEW = [COORquadNEW; 0 0] ;
        
    end
    
    
end