function [Vall,Mintf,MESH,DATA] = FictInterfaceDISP_PLATEquadLIN(PhiDEF,PsiDEFf,MESH,DATAcommon,DATA,Mintf,PhiRB,PsiRBf,DATAoffline)
% Fictitious interface modes for linear quadrilateral plate elements
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/02_PLATES/README_PLATES.mlx
% Partially based on function /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/LARGE_STRAINS/...
% InputDataFunctions/DirichletCONDtime_PLATESfromFEshape.m
% JAHO, 18-FEB-2023
% ----------------------------------------
if nargin == 0
    load('tmp.mat')
    DATAcommon.INCLUDE_FLUCTUATION_DOFS_IN_COARSE_SCALE = 0;
    %    DATAcommon.RESTRICT_NUMBER_FLUCTUATION_MODES_CRITERION = '';
    end

% 1. Coordinate boundary nodes
BoundaryNodes = MESH.faceNODESall ;
%COORbnd = MESH.COOR(BoundaryNodes,:) ;



%LINES_INDEX_FACE = {[MESH.BNDinterface(4).NODES,FACES_BND(1)],[FACES_BND(1),FACES_BND(2)],[FACES_BND(2),FACES_BND(3)],[FACES_BND(3),FACES_BND(4)]} ;

COMB_FACES = {[1,2],[2,3],[3,4],[4,1]};
%COMB_FACES = {[2,3],[3,4],[4,1],[1,2]};

NODES_LINES =  cell(1,4) ;

NODES_CORNERS = cell(size(NODES_LINES)) ;
DeltaZ = zeros(length(NODES_LINES),1) ;

for ilines = 1:length(COMB_FACES)
    iface = COMB_FACES{ilines}(1) ;
    jface = COMB_FACES{ilines}(2) ;
    % THESE ARE THE NODES OF THE EDGES OF THE TRAINING DOMAIN  (LINES)
    NODES_LINES{ilines} =    intersect(MESH.BNDinterface(iface).NODES,MESH.BNDinterface(jface).NODES) ;
    
    % FROM THESE LINES, WE TAKE THE TOP AND BOTTOM NODES (CORNERS). IT IS
    % TACITLY ASSUMED THAT THE TRAINING DOMAIN IS A RECTANGULAR CUBOID,
    % WITH ITS THICKNESS ALONG Z AXIS SIGNIFICANTLY SMALLER THAN IN THE
    % OTHER TWO DIRECTIONS
    [zmax,cornerlocMAX] =  max(MESH.COOR(NODES_LINES{ilines},3)) ;
    [zmin,cornerlocMIN] =  min(MESH.COOR(NODES_LINES{ilines},3)) ;
    DeltaZ(ilines)= zmax-zmin ;
    NODES_CORNERS{ilines} = [NODES_LINES{ilines}(cornerlocMIN),NODES_LINES{ilines}(cornerlocMAX)]' ;
    
end


NODES_CORNERS = cell2mat(NODES_CORNERS) ;
NODES_CORNERS_bottomTOP = NODES_CORNERS;
NODES_CORNERS = NODES_CORNERS' ;
NODES_CORNERS = NODES_CORNERS(:) ;


% FOR LATER PURPOSES, COMPUTATION OF GAUSS POINTS FOR THE HEXAHEDRA ELEMENT
% Patterned after /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/DOMAINdecompHROM/MultiHROM/FictInterfaceDISP_QUADlin.m
% See further comments in /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/02_PLATES/README_PLATES.mlx
% FOR COMPARISON PURPOSES WITH THE CECM INTEGRATION RULE
MESHcoar_COOR = MESH.COOR(NODES_CORNERS,:) ;



MESHLOC.COOR = MESHcoar_COOR  ;
MESHLOC.CN = [1,2,3,4,5,6,7,8];
DATAloc.MESH = MESHLOC ;
DATAloc.NumberOfGaussPointsPerElement = [3,3,3] ;
DATAloc.MESH.TypeElement = 'Hexahedra' ;
[~,DATA.weigthsGAUSSIAN,DATA.xGAUSSIAN] = GetMeshVariablesGAUSS(DATAloc) ;








rnodLOC = unique(BoundaryNodes) ;  % LIST OF NODES WITH PRESCRIBED DISPLACEMENTS (BOUNDARY)
% ----------------------------------------------------------
COORbnd = MESH.COOR(rnodLOC,:) ;  % COORDINATES

% SHAPE FUNCTIONS
% ---------------
% See example in /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/ContinuousEmpiricalCubatureM/Paper_hernandez2021ecm/09_ASSESSMENTpaper/Fun3D/Dpoly3D_P5.m
%

[nnodeBND, ndim] = size(MESHcoar_COOR) ;
% Order of polynomial
nnodes_DIM = (nnodeBND)^(1/ndim)  ;  % Number of nodes per dimension
ORDER_POLYNOMIALS = (nnodes_DIM-1)*ones(1,ndim) ;
DATAshape = ShapeFunCoefficients(MESHcoar_COOR,ORDER_POLYNOMIALS) ;
DATAlocSHAPE.DATAshape  = DATAshape;
xLIM = [] ;
DATAlocSHAPE.ORDER_POLYNOMIALS = ORDER_POLYNOMIALS;
[Nshape, ~,~ ]=    ShapeFunctionFE(xLIM,COORbnd,DATAlocSHAPE) ;

% NEXT STEP: DETERMINE THE "MODES" ASSOCIATED TO THESE SHAPE FUNCTIONS (8x 3 = 24 )
Vfe= zeros(size(Nshape,1)*ndim,size(Nshape,2)*ndim) ;
for idim = 1:ndim
    Vfe(idim:ndim:end,idim:ndim:end) = Nshape  ;
end

% PLATE MODES
Vrb = TransformHEXA_COORD_2_PLATE_COOR(Vfe,DeltaZ) ;



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

 

GidPostProcessModesDOML(COORbnd,CNbREN,MESH.TypeElementB,[Vrb,PhiFLUC],[],NameFileMesh,NameFile_res,[],DATALOC) ;


%%%%%%%%%%% GEOMETRIC INFORMATION MIDPLANE ***********************++++
% COORDINATES MIDPLANE 
COORmid = 0.5*(MESHcoar_COOR(1:4,:) + MESHcoar_COOR(5:8,:)); 
thickness = (MESHcoar_COOR(1:4,3) - MESHcoar_COOR(5:8,3)); 
thickness = sum(abs(thickness))/length(thickness) ; 
FESHAPE_coarse_elem_transf_coord.COOR =COORmid ; % 
FESHAPE_coarse_elem_transf_coord.CN =1:size(COORmid,1) ; 
FESHAPE_coarse_elem_transf_coord.thickness = thickness;

%[nnodeBND, ndim] = size(MESHcoar_COOR) ;
% Order of polynomial
%ORDER_POLYNOMIALS = ones(1,2) ;
%DATAshape = ShapeFunCoefficients(COORmid(:,1:2),ORDER_POLYNOMIALS) ;
%DATAlocSHAPE.DATAshape  = DATAshape;
%xLIM = [] ;
%DATAlocSHAPE.ORDER_POLYNOMIALS = ORDER_POLYNOMIALS;
%[Nshape, ~,~ ]=    ShapeFunctionFE(xLIM,COORbnd,DATAlocSHAPE) ;


FESHAPE_coarse_elem_transf_coord.ShapeFunction =''; %'ShapeFunctionFE' ; 
%FESHAPE_coarse_elem_transf_coord.DATAlocSHAPE =DATAlocSHAPE; 

DATA.FESHAPE_coarse_elem_transf_coord = FESHAPE_coarse_elem_transf_coord; 


% MESHLOC.COOR = MESHcoar_COOR  ;
% MESHLOC.CN = [1,2,3,4,5,6,7,8];
% DATAloc.MESH = MESHLOC ;
% DATAloc.NumberOfGaussPointsPerElement = [3,3,3] ;
% DATAloc.MESH.TypeElement = 'Hexahedra' ;
% [~,DATA.weigthsGAUSSIAN,DATA.xGAUSSIAN] = GetMeshVariablesGAUSS(DATAloc) ;
% 

% FESHAPE_coarse_elem_transf_coord.DATA_ShapeFunctionFE =DATAlocSHAPE ; 
% FESHAPE_coarse_elem_transf_coord.ShapeFunction = 'ShapeFunctionFE' ; 
% 
% 



DATAcommon= DefaultField(DATAcommon,'INCLUDE_FLUCTUATION_DOFS_IN_COARSE_SCALE',0) ;

if DATAcommon.INCLUDE_FLUCTUATION_DOFS_IN_COARSE_SCALE == 0
    Vall = Vrb;
else
    
    
    disp('FILTERING FLUCTUATION MODES ')
    nfluc = size(PhiFLUC,2) ;
    DATAcommon= DefaultField(DATAcommon,'RESTRICT_NUMBER_FLUCTUATION_MODES_CRITERION','BY_ANGLES_Vrb') ;
    
    PhiFLUC = FilterFluctuationModes(DATAcommon,PsiRBf,PhiFLUC,PsiDEFf,Vrb,Mintf) ;
    disp(['Number of included fluctuation modes = ',num2str(size(PhiFLUC,2)),' (of ',num2str(nfluc),')'])
    
    
    % ---------------------------------------------------------------
    % REPLACE THESE MODES BY MODES WITH ZERO VALUE AT THE CORNERS
    % ---------------------------------------------------------------
    PhiFLUC_LINES =  cell(1,4)  ; % FLUCTATION EDGES
    
    ndimPLATE = 5;
    for imodeFLUC = 1:size(PhiFLUC,2)  % Loop over fluctuation modes
        for icorner = 1:size(NODES_CORNERS_bottomTOP,2)  % Loop over corners
            
            % Bottom node,
            cornGLO_bottom = NODES_CORNERS_bottomTOP(1,icorner) ; % Global index corner
            ibot = find(BoundaryNodes == cornGLO_bottom) ; % Local index corner
            dx_bot = PhiFLUC(3*ibot-2,imodeFLUC) ;  %
            dy_bot = PhiFLUC(3*ibot-1,imodeFLUC) ;  %
            dz_bot = PhiFLUC(3*ibot,imodeFLUC) ;  %
            
            % Top node,
            cornGLO_top = NODES_CORNERS_bottomTOP(2,icorner) ; % Global index corner
            itop = find(BoundaryNodes == cornGLO_top) ; % Local index corner
            dx_top= PhiFLUC(3*itop-2,imodeFLUC) ;  %
            dy_top = PhiFLUC(3*itop-1,imodeFLUC) ;  %
            dz_top = PhiFLUC(3*itop,imodeFLUC) ;  %
            
            % PLate displacements
            plate_displ = zeros(1,ndimPLATE) ;
            
            plate_displ(1) = 0.5*(dx_bot+dx_top) ;
            plate_displ(2) = 0.5*(dy_bot+dy_top) ;
            plate_displ(3) = 0.5*(dz_bot+dz_top) ;
            plate_displ(4) =  (dy_bot-dy_top)/DeltaZ(icorner) ;
            plate_displ(5) =  (-dx_bot+dx_top)/DeltaZ(icorner) ; ;
            %
            %
            %
            for idim = 1:ndimPLATE  % Loop over spatial dimension
                imode = ndimPLATE*(icorner-1)+idim ; % Mode corresponding to the corner under consideration
                ValModeAtCorner = plate_displ(idim)  ;
                PhiFLUC(:,imodeFLUC) =  PhiFLUC(:,imodeFLUC)  - ValModeAtCorner*Vrb(:,imode);
            end
            
            
            
            
        end
    end
    
    % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/01_HOMOG/README_HOMOG.mlx
    GidPostProcessModesDOML(COORbnd,CNbREN,MESH.TypeElementB,[Vrb,PhiFLUC],[],NameFileMesh,NameFile_res,[],DATALOC) ;
    
    
    % EACH COLUMN OF PhiFLUC AT THIS STAGE CORRESPONDS TO A FLUCTUATION MODE IN
    % WHICH PLATE DOFS ARE ZERO (NO TRANSLATIONS AND NO ROTATIONS OF THE
    % EDGES). However, compatibility demands that  fluctuations of the
    % edges be considered separately. ACCORDINGLY, NEXT WE DECOMPOSE EACH
    % PhiFLUC(:,imodeFLUC) INTO FIVE DIFFERENT MODES: ONE FOR THE BOUNDARY
    % SURFACES PER SE, AND FOUR MODES FOR THE EDGE FLUCTUATIONS
    % In order not to introduce undesired, irrelevant modes, it is
    % necessary to store edge fluctuations which are large in comparison
    % with the remaining fluctuations
    % -------------------------------------------------------------------------------
    TOL_include_edge_fluctuations = 0.001;  % This is the threshold above which edge fluctuations are incorporated
    % as actual modes
    PhiFLUCedges = cell(1,size(PhiFLUC,2)) ;
    for imodeFLUC = 1:size(PhiFLUC,2)  % Loop over fluctuation modes
        PhiFLUCsurf_loc = PhiFLUC(:,imodeFLUC) ;
        normFLUC_all = max(abs(PhiFLUCsurf_loc)) ;
        PhiFLUCedge_loc = zeros(size(PhiFLUCsurf_loc,1),4) ;
        for icorner = 1:length(NODES_LINES)  % Loop over corners
            [dummy,NodesLOCcorner,~] = intersect(rnodLOC,NODES_LINES{icorner}) ;
            DOFslocLINE = small2large(NodesLOCcorner,ndim) ;
            normFLUCedge = max(abs(PhiFLUCsurf_loc(DOFslocLINE))) ;
            if normFLUCedge/normFLUC_all > TOL_include_edge_fluctuations
                PhiFLUCedge_loc(DOFslocLINE,icorner) =   PhiFLUCsurf_loc(DOFslocLINE) ;
            end
            PhiFLUCsurf_loc(DOFslocLINE) = 0 ;
        end
        PhiFLUC(:,imodeFLUC) = PhiFLUCsurf_loc ;
        TOLlocs = 1e-4 ; DATAlocs.RELATIVE_SVD = 1;
        [Ufluc,Sfluc,Vfluc] = SVDT(PhiFLUCedge_loc,TOLlocs,DATAlocs)  ;
        PhiFLUCedges{imodeFLUC}  = Ufluc ;
        
    end
    
    % FLUCTUATION MODES FOR EDGES
    PhiFLUCedges = cell2mat(PhiFLUCedges) ;
    TOLlocs = 1e-4 ; DATAlocs.RELATIVE_SVD = 1;
    [PhiFLUCedge,Sfluc,Vfluc] = SVDT(PhiFLUCedge_loc,TOLlocs,DATAlocs)  ;
    
    GidPostProcessModesDOML(COORbnd,CNbREN,MESH.TypeElementB,[PhiFLUCedge],[],NameFileMesh,NameFile_res,[],DATALOC) ;
    
    
    
    
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
    GidPostProcessModesDOML(COORbnd,CNbREN,MESH.TypeElementB,[cell2mat(PhiFLUC_faces)],[],NameFileMesh,NameFile_res,[],DATALOC);
    
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
    disp('Master slave strategy (ignoring corners for the time being)')
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
        
        
        GidPostProcessModesDOML(COORbnd,CNbREN,MESH.TypeElementB,[Vall],[],NameFileMesh,NameFile_res,[],DATALOC);
        
        
        
        
    end
    
    
    error('Option not finished yet')
    
end



