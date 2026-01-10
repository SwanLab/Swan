function [Vall,MESH,DATA,INFOPLOTMESHBND] = FictInterfaceDISP_HEXAquad(PhiDEF,PsiDEFf,MESH,DATAcommon,DATA,Mintf,PhiRB,PsiRBf,DATAoffline)
% Fictitious interface modes for hexahedra elements (27 nodes)
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/06_3D_27_HEXAHEDRA/02_27nodeHEX.mlx
% Based on its counterpart for linear elements:
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/01_HOMOG/06_Hexahedra.mlx
% JAHO, 27-MARCH-2023
% ----------------------------------------
if nargin == 0
    load('tmp1.mat')
    
end

% 1. Coordinate boundary nodes
BoundaryNodes = MESH.faceNODESall ;
COORbnd = MESH.COOR(BoundaryNodes,:) ;

%2.) Determination LINES POINTS (as intersections)
PLANES = 1:6 ;
% LINES
LINES = cell(1,12) ;
% % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/06_3D_27_HEXAHEDRA/02_27nodeHEX.mlx
% LINES DEFINED AS INTERSECTION OF PLANES 1,2...6
% RECALL THAT THE ORDER IS (XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX)
LINES{1} = [5,3]   ;
LINES{2} = [5,2]   ;;
LINES{3} = [5,4]   ;
LINES{4} = [5,1]   ;
LINES{5} = [1,3]   ;
LINES{6} = [3,2]   ;
LINES{7} = [2,4]   ;
LINES{8} = [4,1]   ;
LINES{9} =  [6,3]   ;
LINES{10} = [6,2]   ;
LINES{11} = [6,4]   ;
LINES{12} = [6,1]   ;

% POINTS
% -----------------------------------
POINTS.CORNERS = cell(1,8);    % AS INTERSECTIONS OF LINES
%-----------------------------------
POINTS.CORNERS{1} = [4,1] ;
POINTS.CORNERS{2} = [1,2] ;
POINTS.CORNERS{3} = [2,3] ;
POINTS.CORNERS{4} = [3,4] ;
%%----------------------------------
POINTS.CORNERS{5} = [12,9] ;
POINTS.CORNERS{6} = [9,10] ;
POINTS.CORNERS{7} = [10,11] ;
POINTS.CORNERS{8} = [11,12] ;
%-----------------------------------
% POINTS MIDLINES (CENTROID), INDICATE NUMBER OF LINE
% ------------------------------------------------------
POINTS.LINES  = 1:12    ;
% -----------------------

% --------------------------------------------------
% POINTS CENTROID FACES
POINTS.FACES = [5,3,2,4,1,6] ;

%%%%% -----------------------------------------------------------------------------------------


% ------------------------------------------------------------
% NODES OF EACH LINE, as well as centroid (approximated)
% ------------------------------------------------------------
ndim = 3;
COORmidLINES = zeros(length(LINES),ndim) ;
LINES_NODES = cell(size(LINES)) ;
for iline = 1:length(LINES)
    iface = LINES{iline}(1) ;
    jface = LINES{iline}(2) ;
    LINES_NODES{iline} = intersect(MESH.BNDinterface(iface).NODES,MESH.BNDinterface(jface).NODES) ;
    COORLOCnod = MESH.COOR(LINES_NODES{iline} ,:);
    ilineLOC =POINTS.LINES(iline) ;
    COORmidLINES(ilineLOC,:) = sum(COORLOCnod,1)/size(COORLOCnod,1) ;
end
% ------------------------------------------------------------
% CORNER NODES
% ------------------------------------------------------------
COORcorners =zeros(length(POINTS.CORNERS),ndim) ;
for icorner = 1:length(POINTS.CORNERS)
    iline = POINTS.CORNERS{icorner}(1) ;
    jline = POINTS.CORNERS{icorner}(2) ;
    CornerLocInd = intersect(LINES_NODES{iline},LINES_NODES{jline}) ;
    COORcorners(icorner,:) = MESH.COOR(CornerLocInd,:);
end
% ------------------------------------------------------------
% plane NODES  (Centroid, exact)
% ------------------------------------------------------------
COORplanes = zeros(length(POINTS.FACES),ndim) ;
for iplane = 1:length(POINTS.FACES)
    iplaneLOC = POINTS.FACES(iplane) ;
    COORplanes(iplane,:) = MESH.BNDinterface(iplaneLOC).CENTROID ;
end

COORhexa = [COORcorners;COORmidLINES;COORplanes; 0 0 0 ] ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DATAcommon.MESHcoarse
DATAcommon= DefaultField(DATAcommon,'MESHcoarse',[]) ;  % JAHO 18-APR-2024, .MESHcoarse IS A MORE GENERIC WAY TO REFER TO
% THE COARSE, EIFEM ELEMENT. FOR QUADRATIC ELEMENTS IT IS A MUST ! 

DATAoffline= DefaultField(DATAoffline,'MESH_QUADRATIC_PARENT_DOMAIN', []) ;
if  isempty(DATAoffline.MESH_QUADRATIC_PARENT_DOMAIN) && isempty(DATAcommon.MESHcoarse)
    error('yOU MUST DEFINE DATAoffline.MESH_QUADRATIC_PARENT_DOMAIN  OR DATAcommon.MESHcoarse')
end
% Recalculate midside coordinates ---
% 2-Apr-2023: Later versions of the program should use directly the
% information of mesh MESH_QUADRATIC_PARENT_DOMAIN to define the nodes of the coarse-scale element
COORhexa = RecalculatePositionNodesHEXA(COORhexa,DATAoffline.MESH_QUADRATIC_PARENT_DOMAIN,MESH.GEOproperties.CENTROID,DATAcommon.MESHcoarse) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PLOT_POINTS =0;
if PLOT_POINTS == 1
    figure(89)
    hold on
    xlabel('x')
    ylabel('y')
    zlabel('z');
    for ipoints = 1:size(COORhexa)
        plot3(COORhexa(ipoints,1),COORhexa(ipoints,2),COORhexa(ipoints,3),'k*')
        text(COORhexa(ipoints,1),COORhexa(ipoints,2),COORhexa(ipoints,3),['  ',num2str(ipoints)])
    end
end





% FOR COMPARISON PURPOSES WITH THE CECM INTEGRATION RULE
MESHLOC.COOR = COORhexa  ;
MESHLOC.CN = [1:27];
DATAloc.MESH = MESHLOC ;
DATAloc.NumberOfGaussPointsPerElement = [3,3,3] ;
DATAloc.MESH.TypeElement = 'Hexahedra' ;
[~,DATA.weigthsGAUSSIAN,DATA.xGAUSSIAN] = GetMeshVariablesGAUSS(DATAloc) ;

%
% OLD_METHOD = 1;
%
% if OLD_METHOD == 1
%     % Before 18-March-2023
%     nnodeBND = size(COORhexa,1) ;
%
%     nnodes_DIM = (nnodeBND)^(1/ndim)  ;  % Number of nodes per dimension
%     ORDER_POLYNOMIALS = (nnodes_DIM-1)*ones(1,ndim) ;
%     DATAshape = ShapeFunCoefficients(COORhexa,ORDER_POLYNOMIALS) ;
%     DATAlocSHAPE.DATAshape  = DATAshape;
%     xLIM = [] ;
%     DATAlocSHAPE.ORDER_POLYNOMIALS = ORDER_POLYNOMIALS;
%     [Nshape, ~,~ ]=    ShapeFunctionFE(xLIM,COORbnd,DATAlocSHAPE) ;
% else
% Method based on inverse mapping
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/03_HETEROG/03_CURVED.mlx

DATAinp = [] ;
TypeElement = 'Hexahedra' ;

DATAoffline.InterpolationMethod = 'FEnodes'  ;

DATAcommon = DefaultField(DATAcommon,'ORDER_INVERSE_ELEMENT_inverse_mapping',4) ;

DATAoffline = DefaultField(DATAoffline,'ORDER_INVERSE_ELEMENT_inverse_mapping',DATAcommon.ORDER_INVERSE_ELEMENT_inverse_mapping) ;
DATAoffline = DefaultField(DATAoffline,'NPOINTS_ONE_DIRECTION_inverse_mapping',10) ;

%     DATAinp.ORDER_INVERSE_ELEMENT_inverse_mapping = DATAoffline.ORDER_INVERSE_ELEMENT_inverse_mapping;
%     DATAinp.NPOINTS_ONE_DIRECTION_inverse_mapping = DATAoffline.NPOINTS_ONE_DIRECTION_inverse_mapping;  ;

[Nshape,~,~] = ShapeFunDer_inversemapping(COORhexa',COORbnd,TypeElement,DATAoffline) ;
DATAlocSHAPE = []  ;
DATAlocSHAPE.DATAshape  = [];
%end




% These are the 9 shape functions corresponding to the 9 nodes of a
% quadratic element. However, the ninth shape function is the one
% corresponding to the "centroid" node, and therefore, its value at the
% boundary is identically zero ---thus, we have to get rid of the last
% column of Nshape
Nshape = Nshape(:,1:end-1) ;


%
% MESHLOC.COOR = COORhexa  ;
% MESHLOC.CN = [1:8];
% DATAloc.MESH = MESHLOC ;
% DATAloc.NumberOfGaussPointsPerElement = [2,2,2] ;
% DATAloc.MESH.TypeElement = 'Hexahedra' ;
% [~,DATA.weigthsGAUSSIAN,DATA.xGAUSSIAN] = GetMeshVariablesGAUSS(DATAloc) ;
%
%
%
%
% nnodeBND = length(CORNERS) ;
%
% nnodes_DIM = (nnodeBND)^(1/ndim)  ;  % Number of nodes per dimension
% ORDER_POLYNOMIALS = (nnodes_DIM-1)*ones(1,ndim) ;
% DATAshape = ShapeFunCoefficients(COORhexa,ORDER_POLYNOMIALS) ;
% DATAlocSHAPE.DATAshape  = DATAshape;
% xLIM = [] ;
% DATAlocSHAPE.ORDER_POLYNOMIALS = ORDER_POLYNOMIALS;
% [Nshape, ~,~ ]=    ShapeFunctionFE(xLIM,COORbnd,DATAlocSHAPE) ;

% Recall that this EIF element will play the role of "parent element" in
% our coares-scale formulation. Therefore, it is necessary to store the
% shape functions and the derivatives of the shape functions for purposes
% of transforming the coordinates when the physical domain does not
% coincide with the parent domain.

FESHAPE_coarse_elem_transf_coord.COOR =COORhexa ; %
FESHAPE_coarse_elem_transf_coord.CN =1:size(COORhexa,1) ;
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
NameLoc =     'DispIntfhexa' ;
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


% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/01_HOMOG/README_HOMOG.mlx

GidPostProcessModesDOML(COORbnd,CNbREN,MESH.TypeElementB,[Vrb,PhiFLUC],[],NameFileMesh,NameFile_res,[],DATALOC) ;
INFOPLOTMESHBND.COOR = COORbnd;
INFOPLOTMESHBND.CN = CNbREN;
INFOPLOTMESHBND.TypeElement = MESH.TypeElementB;
INFOPLOTMESHBND.posgp = [];
INFOPLOTMESHBND.NameFileMesh = NameFileMesh;



DATAcommon= DefaultField(DATAcommon,'INCLUDE_FLUCTUATION_DOFS_IN_COARSE_SCALE',0) ;

if DATAcommon.INCLUDE_FLUCTUATION_DOFS_IN_COARSE_SCALE == 0  ||  isempty(PhiDEF)
    Vall = Vrb;
else
    
    error('THIS IS UNFINISHED ...')
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


