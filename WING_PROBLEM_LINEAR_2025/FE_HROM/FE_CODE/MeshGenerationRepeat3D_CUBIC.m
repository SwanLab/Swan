function [COOR,CN,TypeElement,CONNECTb,TypeElementB,Materials,MATERIALNEW,DOMAINVAR]= ...
    MeshGenerationRepeat3D_CUBIC(NAMEMSH,DATA,MATERIAL)
% Automatic generation of meshes by tiling copies (translation in one
% direction). See NEW_IDEAS.pdf
% Copy of  MeshGenerationRepeat3D_facesVCURV,  able to deal with variable curvature)
% 14-OCT-202O
% JaHO
%
% OUTPUTS:
% --------
%
% DOMAINVAR.ListElements = ListElementsDom ; % List of elements of each domain
% DOMAINVAR.NODES_faces12 = NODES_faces12 ; % List of connecting nodes of each domain (face 1  and face 2), so that
% the (y,z) coordinates of NODES_faces12{idom,iface} and
% NODES_faces12{jdom,iface} are the same for any "idom" and any "iface"
% ----------------------------------------------------------------
%addpath('DATA_input')
% INPUTS
% ---------

%dbstop('11')
if nargin ==0
    load('tmp2.mat')
    % DATA.TypeUnitCell = 'HEXAHEDRA' ;
end
%%%%
rotGLO =[] ;
nDOMglo = DATA.MakeMeshByRepetition.nDOM ; % number of domains per direction
nDIR = length(nDOMglo) ; % Number of directions
if  ~iscell(DATA.MakeMeshByRepetition.DIRECTION)
    DATA.MakeMeshByRepetition.DIRECTION = { DATA.MakeMeshByRepetition.DIRECTION} ;
end


disp('REpeating cells ...')
nDOM = nDOMglo;
SLICE.NAME =NAMEMSH ;
DATA3D = GeometrySlice(SLICE,DATA) ; % New function, for reading mesh and node faces



DATA = DefaultField(DATA,'angROTATION_FACE',[]) ; % Not really necessary now ....
DATA = DefaultField(DATA,'dX',[]) ;
DATA = DefaultField(DATA,'dY',[]) ;

%
%OOR = DATA3D.COOR ; % Coordinates
CN = DATA3D.CN ; % Connectivities
TypeElement = DATA3D.TypeElement   ;
TypeElementB = DATA3D.TypeElementB ; % Type boundary element
if ~isempty(DATA3D.MaterialType)
    MaterialType = DATA3D.MaterialType ;
else
    MaterialType = ones(size(CN,1),1) ;
end
CONNECTb = DATA3D.CNb ;  % Face connectivities

MATERIALNEW.PLY = MATERIAL.PLY ;  % Set of materials defined in the input data

nmat = length(unique(MaterialType)) ;
if nmat ~=length(MATERIAL.PLY)
    error('The number of materials should coincide with the number of materials defined in the input file')
end
NODESfaces = DATA3D.NODES_FACES ;  % All face nodes (labels from GID)

% 8-JUN-2019. Problems appeared when dealing with multisurface geometries.
% CONNECTb contains element that are not genuine boundary elements, and
% this appears to introduce problems when forming the set of boundary elements
% 28th JUN-2019 --- tHE MODIFICATIO MADE THEN APPEARS TO HAVE RUINED OTHER
% SIMULATIONS.

% WE comment it out for the time being (28-Jun-2019). However, bear in mind
% that it might give rise to problems
% warning('Revise this part ..... ')
% NOT_ELIMINATE = 1 ;
%
% if NOT_ELIMINATE == 0
%     ALLNODESFACES= unique(cell2mat(NODESfaces')) ;  % All nodes included by the user as NODES FACES
%     INCLUDEIT = zeros(size(CONNECTb));
%     for inodes = 1:size(CONNECTb,2)
%         NODESLOC = CONNECTb(:,inodes) ;
%         [~, JJJ,KKK] =  intersect(NODESLOC,ALLNODESFACES) ;
%         INCLUDEIT(JJJ,inodes)  = 1;
%     end
%
%     INCLUDEIT = prod(INCLUDEIT,2)  ;
%     INCLUDEIT = find(INCLUDEIT==1) ;
%     CONNECTb = CONNECTb(INCLUDEIT,:) ;
% end
% We only take into account connectivities of the
% faces specified by NODESfaces (by the user, in GID)
% Face 1 and 2 are those normal to the "repetition" direction (plane YZ)
CONNECTb_faces = cell(1,length(NODESfaces)) ;
for iface = 1:length(NODESfaces)
    [dummy, setBelemLOC]= ElemBnd(CONNECTb,NODESfaces{iface}); % elements face "iface"
    % Connectivities faces f1 and f2
    CONNECTb_faces{iface} = CONNECTb(setBelemLOC,:) ;
    % The above is the connectivity matrix for the nodes of face "iface"
end




% CHECKING THAT NO REPEATED ELEMENTS HAVE BEEN INTRODUCED  (Jan-7-2019)
% --------------------------------------------------------
CONNECTb_CHECk =cell2mat(CONNECTb_faces(:)) ;
[CONNECTbNonRepeated,NODES_repeated ]= CheckCONNECTbNONrepeated(CONNECTb_CHECk) ;
if ~isempty(NODES_repeated)
    warning(' ')
    disp('There is something wrong in the construction of CNb')
    disp('Revise the definition of faces in the GID-preprocess file')
    disp('FACES must be disjoints sets')
    disp(['Check faces corresponding to nodes ',num2str(NODES_repeated')])
    disp(['Generate again the mesh (both interior and boundary) and export the .MSH and .DATA files ' ])
    error(' ')
end

% --------------------------------------------------------


NODES_faces12 = cell(nDOM,2) ;

ndim = length(DATA3D.CENTRf1) ;

if nDOM== 1
    % There is just one single domain
    % --------------------------------
    MATERIALNEW =    MATERIAL   ;
    Materials =  MaterialType   ;
    CONNECTb = cell(nDOM,length(CONNECTb_faces))  ;
    for iface = 1:length(NODESfaces)
        CONNECTb{1,iface} = CONNECTb_faces{iface} ; % Each entry should be a cell itself
    end
    ListElementsDom = {1:size(CN,1)} ;
    ListNodesDom = {1:size(COOR,1)} ;
    NODES_faces12{1,1} = DATA3D.NODES_FACES{1} ;
    NODES_faces12{1,2} = DATA3D.NODES_FACES{2} ;
    
    % Stretching of the domain in the x-direction (with respect centroid face 1)
    DDD = DATA.MakeMeshByRepetition  ;
    DDD = DefaultField(DDD,'DILATATION_SLICES_x',1) ;
    eFACTOR =DDD.DILATATION_SLICES_x(1) ;
    idir = 1;
    DATA3D{1}.COOR(:,idir) = DATA3D{1}.CENTRf1(idir) + eFACTOR*( DATA3D{1}.COOR(:,idir) - DATA3D{1}.CENTRf1(idir)) ;
    
    
else
    
    % Rotation matrix
    % --------------
    %  DATA3D = DefaultField(DATA3D,'rotDOMfacesLOC',[]) ;
    % rotLOC = cell(size(DATA3D)) ;
    % for idom = 1:length(DATA3D)
    %     rotLOC{idom} = DATA3D{idom}.rotDOMfacesLOC{2} ;
    % end
    rotGLO = cell(1,nDOM) ;
    %% iface2 = 2;
    % if ~isempty(DATA3D.rotDOMfacesLOC)
    %    if ~isempty(DATA3D.rotDOMfacesLOC{iface2})
    %       rotLOC = DATA3D.rotDOMfacesLOC{iface2} ;   % Local rotation
    
    if iscell(DATA.ANGiniGLO)
        rotGLO{1}  = DATA.ANGiniGLO{1} ;
    else
        thetaINI  = DATA.ANGiniGLO(1);
        rotGLO{1} = [cos(thetaINI) -sin(thetaINI) 0
            sin(thetaINI) cos(thetaINI) 0
            0              0          1] ;
    end
    %  end
    %end
    
    % -----------
    f1NOD = DATA3D.NODES_FACES{1} ;  % Nodes face 1, domain 1
    f2NOD = DATA3D.NODES_FACES{2} ;  % Nodes face 2, domain 1,
    % Notice that f1NOD(i) matches with f2NOD(i)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Global arrays
    COORglo =DATA.COORtransf;  % Global coordinates, already transformed in  function GeneralCurvature.m
    DATA.COORtransf = [] ;
    
    %     % Translation vector
    %     translationVECTOR= cell(nDOM,1) ;  % Global coordinates
    %
    %     for idom = 1:nDOM
    %         translationVECTOR{idom} =  DATA3D{idom}.CENTRf2- DATA3D{idom}.CENTRf1;  % This should be changed for curved elements
    %     end
    
    
    
    
    CNglo   = cell(nDOM,1) ;  % Global connectivities, all nodes
    CNglo{1} = [CN]  ;
    CONNECTbGLO=  cell(nDOM,length(CONNECTb_faces)) ; % Global connectivities, boundary nodes
    CONNECTbGLO(1,:) =  CONNECTb_faces;
    
    % We create also an array for storing the nodes of faces F1 and F2 of
    % each domain.
    NODES_faces12 = cell(nDOM,2) ;
    NODES_faces12{1,1} = f1NOD ;
    NODES_faces12{1,2} = f2NOD ;
    % ---------------------------
    Materials =cell(nDOM,1) ;
    Materials{1} =MaterialType;
    % --------------------------
    ListElementsDom = cell(nDOM,1) ;
    ListElementsDom{1} = 1:size(CN,1) ;  % List of elements domain 1
    ListNodesDom = cell(nDOM,1) ;
    ListNodesDom{1} = 1:size(DATA3D.COOR,1) ;  % List of nodes domain 1
    
    %% Loop over domains
    [nnodes,ndim] = size(DATA3D.COOR) ;
    
    %%% Coordinate relative to centroid face 1
    %     COOR_REL = COOR ; %
    %     for idim = 1:size(COOR,2)
    %         COOR_REL(:,idim) =    COOR(:,idim) - DATA3D.CENTRf1(idim) ;
    %     end
    rotTOTAL =  rotGLO{1} ;   % total rotation
    for e = 2:nDOM
        disp(['Cell =',num2str(e)])
        
        
        %  -----------------------
        % Coordinates domain e ***
        % ------------------------
        %%% STRETCHING FACTOR  for this domain
        %   eFACTOR = DATA.MakeMeshByRepetition.DILATATION_SLICES_x(e) ;
        % Computation of the global coordinates of this domain
        %  COORglo{e} =  zeros(size(DATA3D{e}.COOR))  ;
        
        
        
        if iscell(DATA.ANGiniGLO)
            rotGLO{e}  = DATA.ANGiniGLO{e} ;
        else
            thetaINI = DATA.ANGiniGLO(e) ;
            
            rotGLO{e} = [cos(thetaINI) -sin(thetaINI) 0
                sin(thetaINI) cos(thetaINI) 0
                0              0          1] ;
        end
        
        
        
        
        
        
        %% List of nodes
        % -----------------------------------------------------
        ListNodesDom{e} = ListNodesDom{e-1}+nnodes ; % List of nodes domain ( for ensuring same numbering
        % of nodes for all domains)
        % Connectivity matrix domain. All elements
        %----------------------------------------------
        CNnew = CNglo{e-1} + nnodes; % We sum up the number of nodes of each domain
        ListElementsDom{e} = ListElementsDom{e-1}+size(CNnew,1) ;
        % Boundary elements (for each labeled face)
        CNbLOC = CONNECTbGLO(e-1,:) ;
        for iface =1:length(CONNECTb_faces)
            CNbLOC{iface} = CONNECTbGLO{e-1,iface} + nnodes;
        end
        % -------------------------------------------
        %         f1NEW = f1NOD + (e-1)*nnodes ;  % Numbering of new face 1
        %         f2OLD = f2NOD + (e-2)*nnodes ;  % Numbering of old face 2
        f1NEW = NODES_faces12{e-1,1} + nnodes;
        f2OLD = NODES_faces12{e-1,2} ;
        % -------------------------------------------------------------
        NODES_faces12{e,1} = f2OLD ;   % Face 1,
        NODES_faces12{e,2} = NODES_faces12{e-1,2} + nnodes ;   % Face 2
        % ----------
        
        % Renumbering: we have to replace nodes (f1NEW) by f2NOD
        % in CNnew and CNbLOC
        ListNodesDom{e}(f1NOD) =  ListNodesDom{e-1}(f2NOD) ;
        
        CNnewREN  = CNnew ;
        CNbLOCnew = CNbLOC ;
        for ifacen = 1:length(f1NOD)   % This may be improved (make it more efficient)...
            nodeLOC = f1NEW(ifacen);
            % Replacing it in the connectivity matrix
            INDnodes = find(CNnew==nodeLOC) ;
            CNnewREN(INDnodes) = f2OLD(ifacen) ;
            % Replacing it in the face connectivity matrix
            for iface = 1:length(CNbLOCnew)
                INDnodesb = find(CNbLOC{iface}==nodeLOC) ;
                CNbLOCnew{iface}(INDnodesb) =f2OLD(ifacen) ;
            end
        end
        
        CNglo{e} =  CNnewREN ;
        CONNECTbGLO(e,:) = CNbLOCnew ; %
        Materials{e} = Materials{e-1} + nmat ;
        for imat = 1:length(MATERIAL.PLY)
            matREF = (e-1)*nmat ;
            MATERIALNEW.PLY(matREF+imat) =  MATERIAL.PLY(imat) ;
        end
        
    end
    
    Materials = cell2mat(Materials) ;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('REnumbering...')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Numeration of nodes is not consecutive. We turn it consecutive in
    % what follows
    CNglo = cell2mat(CNglo) ;
    NODES =  unique(CNglo);  % list of "old" node numbers, sorted from small to large
    COORglo= cell2mat(COORglo') ;
    
    COOR = COORglo(NODES,:) ;          % Global matrix of coordinates
    NODES_new = 1:length(NODES) ;      % List of "new" node numbers, from 1 to length(NODES)
    % Interior CN
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CN = RenumberConnectivities(CNglo,NODES_new) ;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Boundary CN
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CONNECTb = CONNECTbGLO ;
    for idom = 1:size(CONNECTbGLO,1)
        for iface = 1:size(CONNECTbGLO,2)
            NODESbnd = unique(CONNECTbGLO{idom,iface}) ;
            [~,NODES_bnd,~] = intersect(NODES,NODESbnd) ;
            CONNECTb{idom,iface}= RenumberConnectivities(CONNECTbGLO{idom,iface},NODES_bnd) ;
        end
        
    end
    
    %% List of nodes FACES 1 and 2
    for idom = 1:nDOM
        for iface =1:size(NODES_faces12,2)
            [~,NODES_bnd,~] = intersect(NODES,NODES_faces12{idom,iface}) ;
            NODES_faces12{idom,iface}= RenumberConnectivities(NODES_faces12{idom,iface},NODES_bnd) ;
        end
    end
    
    %% List of all nodes
    for idom = 1:nDOM
        [~,NODES_dom,~] = intersect(NODES,ListNodesDom{idom}) ;
        ListNodesDom{idom}= RenumberConnectivities(ListNodesDom{idom},NODES_dom) ;
    end
    
end

DOMAINVAR.ListElements = ListElementsDom ; % List of elements of each domain
DOMAINVAR.NODES_faces12 = NODES_faces12 ; % List of connecting nodes of each domain (face 1  and face 2), so that
% the (y,z) coordinates of NODES_faces12{idom,iface} and
% NODES_faces12{jdom,iface} are the same for any "idom" and any "iface"
DOMAINVAR.ListNodesDom = ListNodesDom ;
DOMAINVAR.rotGLO = rotGLO ;




DATA = DefaultField(DATA,'PREPRINT_MESH',1) ;

if DATA.PREPRINT_MESH == 1
    d = [] ; strainGLOgid = [] ;  stressGLOgid = []; React = [] ; posgp = [] ;
    NameFileMesh = ['prueba.msh'] ;
    Fnodes = [] ;
    
    NAME_BASE_GIDfiles = DATA.INPUTDATAfile ; %fileparts() ; 
    
    DATA.MATERIAL = MATERIAL ; 
    
    GidPostProcess(COOR,CN,TypeElement,d,strainGLOgid, stressGLOgid,  ...
        React,NAME_BASE_GIDfiles,posgp,NameFileMesh,Materials,DATA,Fnodes);
     
    error('Disable this error command to go on !!!!')
    
end




end

