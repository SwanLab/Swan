function COORhexaNEW = RecalculatePositionNodesHEXA20(COORhexa,MESH_QUADRATIC_PARENT_DOMAIN,CENTROID) ;
% Recalculation of midside, midplane nodal coordinates, quadratic hexahedra  
% element
% Patterned after MidsideCoordintesQQ.m 
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/06_3D_27_HEXAHEDRA/02_27nodeHEX.mlx
% JAHO, 27-March-2023
%------------------------------
if nargin == 0
    load('tmp.mat')
end
COORhexaNEW = COORhexa;
if ~isempty(MESH_QUADRATIC_PARENT_DOMAIN)
    
    
    %% THESE ARE THE COORDINATES OF THE POINTS READ FROM THE TRAINING MESH FILE
    COORgidPARENT = COORhexa ;
    ndim = size(COORgidPARENT,2);
    
    for idim = 1:ndim
        COORgidPARENT(:,idim) = COORgidPARENT(:,idim) + CENTROID(idim) ;
    end
    % LET US REFER THEM TO ITS PSEUDO-CENTROID (computed using the corners, not the midsides)
   % COORgidPARENT  = COORgidPARENT(1:end-1,:) ;
    CENTparent = sum(COORgidPARENT(1:8,:),1)/8;
    for idim = 1:ndim
        COORgidPARENT(:,idim) = COORgidPARENT(:,idim) - CENTparent(idim) ;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % NOW WE READ THE COORDINATES OF THE PARENT DOMAIN CREATED BY GID
    MESHparent= ReadMeshFileStr(MESH_QUADRATIC_PARENT_DOMAIN,'READ_MATERIAL_COLUMN',0)  ;
    
    COORgid = MESHparent.COOR(MESHparent.CN,:) ;
%     if size(MESHparent.COOR,1) == 27
%         COORgid = COORgid(1:end-1,:) ;
%     end
    % LET US REFER THEM TO ITS PSEUDO-CENTROID
    CENT = sum(COORgid(1:8,:),1)/8;
    for idim = 1:ndim
        COORgid(:,idim) = COORgid(:,idim) - CENT(idim) ;
    end
    
    k = dsearchn(COORgid,COORgidPARENT) ;
    
    DISTANCE= COORgidPARENT -COORgid(k,:)  ;
    
    
    distCORNERSnorm = norm(DISTANCE(1:8,:)) ;
    
         dtyp = norm(COORgidPARENT(:,1)-COORgidPARENT(:,2)) ; 
     TOL = 1e-6; 
%     if distCORNERSnorm > TOL*dtyp
%         error('Incorrect parent domain')
%     else

    
    if distCORNERSnorm >  TOL*dtyp
        error('Incorrect parent domain')
    else
        % Therefore
        COORhexaNEW = COORgid(k,:)  ; ;
        for idim = 1:ndim
            COORhexaNEW(:,idim) = COORhexaNEW(:,idim) + CENTparent(idim) - CENTROID(idim) ;
        end
      %  COORhexaNEW = [COORhexaNEW; 0 0 0] ;
        
    end
    
    
end