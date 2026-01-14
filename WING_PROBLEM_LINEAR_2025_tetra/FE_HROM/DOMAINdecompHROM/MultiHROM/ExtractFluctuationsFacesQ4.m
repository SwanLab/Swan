function  FluctuationFacesGLO =  ExtractFluctuationsFacesQ4(Vrb,PhiDEF,MESH,DATAcommon,DATA,Mintf,PhiRB,PsiRBf,DATAoffline,...
    CORNERS,CNbREN,NameFileMesh,NameFile_res,DATALOC,PsiSEf)

if nargin == 0
    load('tmp.mat')
end

disp(' FLUCTUATION MODES (FACE-WISE ) ')
BoundaryNodes = MESH.faceNODESall ;
COORbnd = MESH.COOR(BoundaryNodes,:) ;

PhiDEF_f = PhiDEF(MESH.faceDOFSall,:) ;
nfaces = 4; 
 FluctuationFaces = cell(1,nfaces) ;  
 ndim  = size(MESH.COOR,2) ; 
 
FluctationsTOGETHER =PhiDEF_f ; 

for imodeFLUC = 1:size(PhiDEF_f,2)
    for innode = 1:length(CORNERS)
        cornGLO = CORNERS(innode) ; % Global index corner
        cornLOC = find(cornGLO == BoundaryNodes) ; % Local index corner
        %   dofCORNloc = small2large(cornLOC,ndim) ;
        for idim = 1:ndim  % Loop over spatial dimension
            imode = ndim*(innode-1)+idim ; % Mode corresponding to the corner under consideration
            idof = ndim*(cornLOC-1) +idim ; % DOF corresponding to the corner/dimension under consideration
            ValModeAtCorner = PhiDEF_f(idof,imodeFLUC)  ;
            FluctationsTOGETHER(:,imodeFLUC) =  FluctationsTOGETHER(:,imodeFLUC)  - ValModeAtCorner*Vrb(:,imode);
        end
    end
end

% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/01_HOMOG/README_HOMOG.mlx
GidPostProcessModesDOML(COORbnd,CNbREN,MESH.TypeElementB, [Vrb,FluctationsTOGETHER],DATA.MESH.posgp,NameFileMesh,NameFile_res,[],DATALOC) ;


disp('WORK DONE BY THE FLUCTUATIONS')
PsiSEf'*FluctationsTOGETHER

% -----------------------------------------------------------------------
%  IN ORDER TO PROPERLY IMPOSE COMPATIBILITY, IT IS NECESSARY TO CONSIDER
%  THESE MODES "SEPARATELY", FACE BY FACE
% -----------------------------------------------------------------------
nmodesFLUC = size(FluctationsTOGETHER,2)         ;
nfaces  =  length(MESH.BNDinterface) ;
FluctuationFaces = cell(1,nfaces)       ;

FluctuationFacesGLO = cell(1,nfaces)       ;


for iface  = 1:nfaces
    FluctuationFaces{iface} = zeros(size(FluctationsTOGETHER,1),nmodesFLUC) ;
    FluctuationFacesGLO{iface} = sparse(size(PhiDEF,1),nmodesFLUC) ;
    
    NodesFaceGlo = MESH.BNDinterface(iface).NODES;  % Global numbering nodes of this face
    DOFsGLO =  small2large(NodesFaceGlo,ndim) ;
    [~,NodesFaceLoc,~] = intersect(MESH.faceNODESall,NodesFaceGlo) ;
    MESH.BNDinterface(iface).NODESloc = NodesFaceLoc ;
    DOFsLOC = small2large(NodesFaceLoc,ndim) ;
    MESH.BNDinterface(iface).DOFsLoc = DOFsLOC ;
    
    FluctuationFaces{iface}(DOFsLOC,:) =  FluctationsTOGETHER(DOFsLOC,:);
    FluctuationFacesGLO{iface}(DOFsGLO,:) =  FluctationsTOGETHER(DOFsLOC,:);
    
    %     [Uu,Su,Vu] = SVDT( PhiDEF_f(DOFsLOC,:))
    
end

GidPostProcessModesDOML(COORbnd,CNbREN,MESH.TypeElementB,cell2mat(FluctuationFaces),DATA.MESH.posgp,NameFileMesh,NameFile_res,[],DATALOC);

% Now we have to investigate whether the modes of opposite faces "match"
% each other ---concept of intersection between subspaces

FACES_match = {[1,3],[2,4]}  ;

for imatch = 1:length(FACES_match)
    iface = FACES_match{imatch}(1) ;
    COOR_i = MESH.BNDinterface(iface).COORrelA_global ;  % Coordinates relative to the centroids of each face
    DOFS_i = MESH.BNDinterface(iface).DOFsLoc ;
    Mintf_i_scal = MESH.BNDinterface(iface).GeometricMassMatrix   ;
    Mintf_i = zeros(ndim*size(Mintf_i_scal,1)) ;
    
    for idim = 1:ndim
        Mintf_i(idim:ndim:end,idim:ndim:end) = Mintf_i_scal ;   % Geometric mass matrix interface
    end
    
    jface = FACES_match{imatch}(2) ;
    COOR_j = MESH.BNDinterface(jface).COORrelA_global ;
    DOFS_j = MESH.BNDinterface(jface).DOFsLoc ;
    
    Idx_j = knnsearch(COOR_j,COOR_i) ;
    
    COOR_j_i = COOR_j(Idx_j,:) ;
    
    errorCOOR = norm(COOR_j_i-COOR_j,'fro') ;
    disp('---------------------------------------------')
    disp(['Checking matching interfaces = ',num2str(iface),' and ',num2str(jface)]) ;
    disp(['ERROR_match = ',num2str(errorCOOR)])
    disp('---------------------------------------------')
    
    Idx_DOFs_j_i = small2large(Idx_j,ndim) ;
    DOFS_j_i = DOFS_j(Idx_DOFs_j_i) ;
    
    
    PhiDEF_f_i = FluctuationFaces{iface}(DOFS_i,:);
    PhiDEF_f_j = FluctuationFaces{jface}(DOFS_j_i,:);
    
    TOLcosINTfluc = 1e-3 ;
    % WdefCAND =   WSVDT([PhiDEF_f{1},PhiDEF_f{2}],Mintf{1},DATALOC) ;
    
    %  Mintf_i = eye(size(Mintf_i)) ;
    [COSINE_ANGLES,WdefCAND]= PRANGLES(PhiDEF_f_i,PhiDEF_f_j,Mintf_i,1-TOLcosINTfluc,1e20) ;
    disp(['COSINES PRINCIPLE ANGLES'])
    COSINE_ANGLES
%     
%     % COMPUTE INTERSECTION
%     
%     PhiDEF_f_faces{iface} = zeros(size(PhiDEF_f_faces{iface},1),size(WdefCAND,2)) ;
%     PhiDEF_f_faces{iface}(DOFS_i,:)  = WdefCAND ;
%     
%     PhiDEF_f_faces{jface} = zeros(size(PhiDEF_f_faces{jface},1),size(WdefCAND,2)) ;
%     PhiDEF_f_faces{jface}(DOFS_j_i,:)  = WdefCAND ;
%     
%     
    
    
    
    
    
end

% plot again subspaces
% 
% GidPostProcessModesDOML(COORbnd,CNbREN,MESH.TypeElementB,[Vrb,cell2mat(PhiDEF_f_faces)],DATA.MESH.posgp,NameFileMesh,NameFile_res,[],DATALOC);
% 
%  