function [Vall,Mintf,MESH,DATA,FluctuationFacesGLO] = FictInterfaceDISP_QUADlin_ONLY_GID(PhiDEF,PsiDEFf,MESH,DATAcommon,DATA,Mintf,PhiRB,PsiRBf,DATAoffline)

if nargin == 0
    load('tmp.mat')
    DATAcommon.INCLUDE_FLUCTUATION_DOFS_IN_COARSE_SCALE = 1;
end

% 1. Coordinate boundary nodes
BoundaryNodes = MESH.faceNODESall ;
COORbnd = MESH.COOR(BoundaryNodes,:) ;
DATAcommon = DefaultField(DATAcommon,'MESH_PARENT_DOMAIN_COARSE',[]) ; 
if isempty(DATAcommon.MESH_PARENT_DOMAIN_COARSE)
    error('yOU MUST DEFINE   DATAcommon.MESH_PARENT_DOMAIN_COARSE (MESH OF THE FINITE ELEMENT ENCLOSING THE SUBDOMAIN UNDER STUDY)')
end

 [MESHparent]= ReadMeshFileStr(DATAcommon.MESH_PARENT_DOMAIN_COARSE ,'READ_MATERIAL_COLUMN',0)  ; 
 
COORlin = MESHparent.COOR(MESHparent.CN,:) ;
CENTROID = MESH.GEOproperties.CENTROID ;
ndim = size(CENTROID,2);
for idim = 1:ndim
    COORlin(:,idim) = COORlin(:,idim) - CENTROID(idim) ;
end

     


nnodeBND = 4;

nnodes_DIM = (nnodeBND)^(1/ndim)  ;  % Number of nodes per dimension
ORDER_POLYNOMIALS = (nnodes_DIM-1)*ones(1,ndim) ;
DATAshape = ShapeFunCoefficients(COORlin,ORDER_POLYNOMIALS) ;
DATAlocSHAPE.DATAshape  = DATAshape;
xLIM = [] ;
DATAlocSHAPE.ORDER_POLYNOMIALS = ORDER_POLYNOMIALS;
[Nshape, ~,~ ]=    ShapeFunctionFE(xLIM,COORbnd,DATAlocSHAPE) ;

% Recall that this EIF element will play the role of "parent element" in
% our coares-scale formulation. Therefore, it is necessary to store the
% shape functions and the derivatives of the shape functions for purposes
% of transforming the coordinates when the physical domain does not
% coincide with the parent domain.

FESHAPE_coarse_elem_transf_coord.COOR =COORlin ; %
FESHAPE_coarse_elem_transf_coord.CN =1:size(COORlin,1) ;
FESHAPE_coarse_elem_transf_coord.DATA_ShapeFunctionFE =DATAlocSHAPE ;
FESHAPE_coarse_elem_transf_coord.ShapeFunction = 'ShapeFunctionFE' ;


DATA.FESHAPE_coarse_elem_transf_coord = FESHAPE_coarse_elem_transf_coord;





nmodes = ndim*size(Nshape,2) ;
ndofsLOC = ndim*size(Nshape,1) ;
Vrb =zeros(ndofsLOC,nmodes) ;

for innode = 1:size(Nshape,2)
    for idim = 1:ndim
        imode = ndim*(innode-1)+idim ;
        Vrb(idim:ndim:end,imode) = Nshape(:,innode) ;
    end
end


% PLOT  MODES
% -----------------
% Coordinates
% COORbnd = MESH.COOR(BoundaryNodes,:) ;
% Connectivities
CNb = cell(length(MESH.BNDinterface),1) ;
for iface  = 1:length(MESH.BNDinterface)
    CNb{iface} = MESH.BNDinterface(iface).CNb ;
end
CNb = cell2mat(CNb(:)) ;
CNbREN  =  RenumberConnectivities( CNb,1:length(BoundaryNodes) );


% ------------------------------------------------------------------------------------------------------------------------
NameLoc =     'DispIntfQuad' ;
NAME_MODES_FOLDER = [cd,filesep,'MODES',filesep];
NAME_MODES_DISP = [NAME_MODES_FOLDER,DATA.NAME_BASE,NameLoc] ;
NameFileMesh = [NAME_MODES_DISP,'.msh'] ;
NameFile_res = [NAME_MODES_DISP,'.res'] ;
DATALOC = [];
% ------------------------------------------

% -----------------
% FLUCTUATIONS
% ------------------
% ORTHOGONAL TO Vrb
if ~isempty(PhiDEF)
    PhiDEF_f = PhiDEF(MESH.faceDOFSall,:) ;
    PhiFLUC =   PprojDEF_operator(Vrb,Mintf,PhiDEF_f);
else
    PhiFLUC = [] ;
 end
% 
% PhiDEF_f = PhiDEF(MESH.faceDOFSall,:) ;
% PhiFLUC =   PprojDEF_operator(Vrb,Mintf,PhiDEF_f);

% ---------------------------------------------------------------
% REPLACE THESE MODES BY MODES WITH ZERO VALUE AT THE CORNERS
% ---------------------------------------------------------------
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/01_HOMOG/README_HOMOG.mlx
GidPostProcessModesDOML(COORbnd,CNbREN,MESH.TypeElementB,[Vrb,PhiFLUC],DATA.MESH.posgp,NameFileMesh,NameFile_res,[],DATALOC) ;


% EXTRACTING FLUCTUATIONS FACE-WISE  (FROM PhiDEF)
EXTRACT_Fluct = 0; 
if EXTRACT_Fluct == 1
 FluctuationFacesGLO =  ExtractFluctuationsFacesQ4(Vrb,PhiDEF,MESH,DATAcommon,DATA,Mintf,PhiRB,PsiRBf,DATAoffline,CORNERS,CNbREN,...
     NameFileMesh,NameFile_res,DATALOC,PsiDEFf) ; 
else
    FluctuationFacesGLO = [] ; 
end





DATAcommon= DefaultField(DATAcommon,'INCLUDE_FLUCTUATION_DOFS_IN_COARSE_SCALE',0) ;

if DATAcommon.INCLUDE_FLUCTUATION_DOFS_IN_COARSE_SCALE == 0
    Vall = Vrb;
else
    
    
    disp('FILTERING FLUCTUATION MODES ')
    nfluc = size(PhiFLUC,2) ;
    DATAcommon= DefaultField(DATAcommon,'RESTRICT_NUMBER_FLUCTUATION_MODES_CRITERION','') ;
    
    PhiFLUC = FilterFluctuationModes(DATAcommon,PsiRBf,PhiFLUC,PsiDEFf,Vrb,Mintf) ;
    disp(['Number of included fluctuation modes = ',num2str(size(PhiFLUC,2)),' (of ',num2str(nfluc),')'])
    
    
    for imodeFLUC = 1:size(PhiFLUC,2)
        for innode = 1:length(CORNERS)
            cornGLO = CORNERS(innode) ; % Global index corner
            cornLOC = find(cornGLO == BoundaryNodes) ; % Local index corner
            %   dofCORNloc = small2large(cornLOC,ndim) ;
            for idim = 1:ndim  % Loop over spatial dimension
                imode = ndim*(innode-1)+idim ; % Mode corresponding to the corner under consideration
                idof = ndim*(cornLOC-1) +idim ; % DOF corresponding to the corner/dimension under consideration
                ValModeAtCorner = PhiFLUC(idof,imodeFLUC)  ;
                PhiFLUC(:,imodeFLUC) =  PhiFLUC(:,imodeFLUC)  - ValModeAtCorner*Vrb(:,imode);
            end
        end
    end
    
    % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/01_HOMOG/README_HOMOG.mlx
    GidPostProcessModesDOML(COORbnd,CNbREN,MESH.TypeElementB,[Vrb,PhiFLUC],DATA.MESH.posgp,NameFileMesh,NameFile_res,[],DATALOC) ;
    
    % -----------------------------------------------------------------------
    %  IN ORDER TO PROPERLY IMPOSE COMPATIBILITY, IT IS NECESSARY TO CONSIDER
    %  THESE MODES "SEPARATELY", FACE BY FACE
    % -----------------------------------------------------------------------
    nmodesFLUC = size(PhiFLUC,2)         ;
    nfaces  =  length(MESH.BNDinterface) ;
    PhiFLUC_faces = cell(1,nfaces)       ;
    
    
    for iface  = 1:nfaces
        PhiFLUC_faces{iface} = zeros(size(PhiFLUC,1),nmodesFLUC) ;
        NodesFaceGlo = MESH.BNDinterface(iface).NODES;  % Global numbering nodes of this face
        [~,NodesFaceLoc,~] = intersect(MESH.faceNODESall,NodesFaceGlo) ;
        MESH.BNDinterface(iface).NODESloc = NodesFaceLoc ;
        DOFsLOC = small2large(NodesFaceLoc,ndim) ;
        MESH.BNDinterface(iface).DOFsLoc = DOFsLOC ;
        
        PhiFLUC_faces{iface}(DOFsLOC,:) =  PhiFLUC(DOFsLOC,:);
        
        %     [Uu,Su,Vu] = SVDT( PhiFLUC(DOFsLOC,:))
        
    end
    
    % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/01_HOMOG/README_HOMOG.mlx
    GidPostProcessModesDOML(COORbnd,CNbREN,MESH.TypeElementB,[Vrb,cell2mat(PhiFLUC_faces)],DATA.MESH.posgp,NameFileMesh,NameFile_res,[],DATALOC);
    
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
        
        
        PhiFLUC_i = PhiFLUC_faces{iface}(DOFS_i,:);
        PhiFLUC_j = PhiFLUC_faces{jface}(DOFS_j_i,:);
        
        TOLcosINTfluc = 1e-3 ;
        % WdefCAND =   WSVDT([PhiFLUC{1},PhiFLUC{2}],Mintf{1},DATALOC) ;
        
        %  Mintf_i = eye(size(Mintf_i)) ;
        [COSINE_ANGLES,WdefCAND]= PRANGLES(PhiFLUC_i,PhiFLUC_j,Mintf_i,1-TOLcosINTfluc,1e20) ;
        disp(['COSINES PRINCIPLE ANGLES'])
        COSINE_ANGLES
        
        % COMPUTE INTERSECTION
        
        PhiFLUC_faces{iface} = zeros(size(PhiFLUC_faces{iface},1),size(WdefCAND,2)) ;
        PhiFLUC_faces{iface}(DOFS_i,:)  = WdefCAND ;
        
        PhiFLUC_faces{jface} = zeros(size(PhiFLUC_faces{jface},1),size(WdefCAND,2)) ;
        PhiFLUC_faces{jface}(DOFS_j_i,:)  = WdefCAND ;
        
        
        
        
        
        
        
    end
    
    % plot again subspaces
    
    GidPostProcessModesDOML(COORbnd,CNbREN,MESH.TypeElementB,[Vrb,cell2mat(PhiFLUC_faces)],DATA.MESH.posgp,NameFileMesh,NameFile_res,[],DATALOC);
    
    %% WORK DONE BY THE basic  MODES
    % --------------------------------------
    Wcand = cell2mat(PhiFLUC_faces);
    WorkIntf = PsiDEFf'*Wcand ;
    
    %%% MASTER SLAVE STRATEGY
    numberDOMAINmodes = size(PhiDEF,2) + size(PhiRB,2) ;
    numberRBmodes = size(Vrb,2) ;
    
    if numberDOMAINmodes ~= numberRBmodes
        error('Option not implemented yet')
        %
    else
        % Master-slave strategy
        % See  /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/01_HOMOG/README_HOMOG.mlx
        % INGREDIENTS:
        Phi = [PhiRB(MESH.faceDOFSall,:),PhiDEF(MESH.faceDOFSall,:)] ;
        % Mass matrix: Mintf
        % All interface modes available
        W =  [Vrb,Wcand] ;
        % Minimization problem ---> Find alpha such that  alpha = min arg_{alpha} ||   Phi*q -W*alpha  ||_Mintf
        % Solution  alpha = (W'*Mintf*W)^{-1}*W^T*Mintf*Phi = Y*q
        Y = (W'*Mintf*W)\(W'*Mintf*Phi) ;
        % Now we have  a = Y q
        % Let us split Y as follows  --->  [a_M; a_S]^T  = [Y(master,:) ; Y(slave,:)] q
        
        
        indMASTER = 1:numberRBmodes ;
        indSLAVE = numberRBmodes+1:size(Y,1) ;
        Ym = Y(indMASTER,:) ;
        Ys = Y(indSLAVE,:) ;
        
        %
        J_sm = Ys*pinv(Ym) ;
        % Therefore
        % u = Vrb*aRB + Vfluc*aFLUC = Vrb*aRB + Vfluc*J_sm*aRB = (Vrb + Vfluc*J_sm) ;
        Vall = Vrb+ Wcand*J_sm ;
        
        
        GidPostProcessModesDOML(COORbnd,CNbREN,MESH.TypeElementB,[Vall],DATA.MESH.posgp,NameFileMesh,NameFile_res,[],DATALOC);
        
        %  DATAcommon= DefaultField(DATAcommon,'INCLUDE_FLUCTUATION_DOFS_IN_COARSE_SCALE',0) ;
        
        %     if DATAcommon.INCLUDE_FLUCTUATION_DOFS_IN_COARSE_SCALE == 0
        %         Vall = Vrb;
        %     else
        %         error('Option not implemented yet')
        %
        %     end
        
        
        
    end
    
end


