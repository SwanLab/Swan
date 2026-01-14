function [COOR,CN,TypeElement,CONNECTb,TypeElementB,Materials,MATERIALNEW,DOMAINVAR,DATA]= ...
    MeshGenerationRepeat3D_CUBICDIFF(NAMEMSH,DATA,MATERIAL)
% See MeshGenerationRepeat3D_CUBICDIFF.xoj/pdf
% Automatic generation of meshes by tiling copies (translation in one
% direction). See NEW_IDEAS.pdf
% Copy of  MeshGenerationRepeat3D_CUBICDIFF,  able to deal with variable curvature,
% and slices with different reference meshes)
% Started 1-NOV-202O, plane from Sofia to Barcelona 
% JaHO
%------------------------------------------------------------------------------------
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
% Loop over reference domains
nrefDOM = length(NAMEMSH) ; % Number of reference domains
DATA3D = cell(1,nrefDOM) ;
f1NOD = cell(1,nrefDOM) ;
f2NOD = cell(1,nrefDOM) ;

CONNECTb_faces = cell(1,nrefDOM) ;
irefFACE = 1; % REference face --- face 1
TYPEslices = DATA.TYPEslices ; 
nmat = zeros(nrefDOM) ; 

% OPERATION ON THE SET OF REFERENCE DOMAINS 
% -----------------------------------------
NFACES_REF= zeros(1,nrefDOM) ; 
for irefdom = 1:nrefDOM
    SLICE.NAME =NAMEMSH{irefdom}; % Name of the mesh of reference domain "irefdom"
    DATA3D{irefdom}  = GeometrySlice(SLICE,DATA) ; % Reading mesh and node faces
   % MATERIALNEW(irefdom).PLY = MATERIAL(irefdom).PLY ;  % Set of materials defined in the input data    
    nmat(irefdom) = length(unique(DATA3D{irefdom}.MaterialType)) ;
    if nmat(irefdom) ~=length(MATERIAL(irefdom).PLY)
        error('The number of materials should coincide with the number of materials defined in the input file')
    end   
     % DATA3D{irefdom}.NODES_FACES contains the nodes of the faces of each
    % reference domain . See *.xoj document 
    
    if irefdom == 1
        TOL = 1e-1 ; % Local tolerance
        CNb1 = DATA3D{irefdom}.CNb(1,1:2) ;
        COOR1  = DATA3D{irefdom}.COOR(CNb1,:) ;
        dCOOR = norm(COOR1(1,:)-COOR1(2,:)) ;
        TOL = TOL*dCOOR ; % Abs. tolerance
        COORrefFACE1 = DATA3D{irefdom}.COOR(DATA3D{irefdom}.NODES_FACES{irefFACE},:) ;
        Cf1 = DATA3D{irefdom}.CENTRf1 ;
        for iii= 1:length(Cf1)
            COORrefFACE1(:,iii) =  COORrefFACE1(:,iii)  - Cf1(iii) ; 
        end
    else
        % Sorting nodes of face 1 of all the reference domains
        COORface1LOC = DATA3D{irefdom}.COOR(DATA3D{irefdom}.NODES_FACES{1},:) ;
        Cf1 = DATA3D{irefdom}.CENTRf1 ;
        for iii= 1:length(Cf1)
            COORface1LOC(:,iii) =  COORface1LOC(:,iii)  - Cf1(iii) ; 
        end
        [IDX DISTANCES]= knnsearch(COORface1LOC,COORrefFACE1) ;
        if any(find(DISTANCES > TOL))
            error('Non-conforming meshes between reference domains . Check tolerance TOL')
        end
         DATA3D{irefdom}.NODES_FACES{1} =  DATA3D{irefdom}.NODES_FACES{1}(IDX) ;
         DATA3D{irefdom}.NODES_FACES{2} =  DATA3D{irefdom}.NODES_FACES{2}(IDX) ;
    end
    
      f1NOD{irefdom} = DATA3D{irefdom}.NODES_FACES{1} ;  % Nodes face 1, 
    f2NOD{irefdom} = DATA3D{irefdom}.NODES_FACES{2} ;  % Nodes face 2,  
    
    CONNECTb_faces{irefdom} = cell(1,length(DATA3D{irefdom}.NODES_FACES)) ;
    % DATA3D{irefdom}.NODES_FACES contains the nodes of the faces of each
    % reference domain         
    NFACES_REF(irefdom) =  length(DATA3D{irefdom}.NODES_FACES) ;  ; % Number of faces per domain (max)
    for iface = 1:length(DATA3D{irefdom}.NODES_FACES)  
        [dummy, setBelemLOC]= ElemBnd(DATA3D{irefdom}.CNb,DATA3D{irefdom}.NODES_FACES{iface}); % elements face "iface"
        CONNECTb_faces{irefdom}{iface} =DATA3D{irefdom}.CNb(setBelemLOC,:) ;
        % Connectivity matrix of each set of boundary nodes (face1, face2... and so on )
    end        
end


TypeElement = DATA3D{1}.TypeElement   ;  % We assume it is the same for all elements...
TypeElementB = DATA3D{1}.TypeElementB ; % Type boundary element
nDOM = length(DATA.COORtransf) ; 

NODES_faces12 = cell(nDOM,2) ;
ndim = length(DATA3D{1}.CENTRf1) ;
if nDOM== 1    
    error('Use a standadard mesh construction routine for this problem ')
else       
    rotGLO = cell(1,nDOM) ;  % Rotation matrix R_G1 of each slice (G = Global Ref. system. 1= Local ref. system at face 1)
    rotGLO{1}  = DATA.ANGiniGLO{1} ;   % Genererated by ParamCurv_DIFFSLICES_3D.m 
    % -----------  
    % Notice that f1NOD(i) matches with f2NOD(i)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Global arrays
    COORglo =DATA.COORtransf;  % Global coordinates, already transformed in  function GeneralCurvature.m
    DATA.COORtransf = [] ;   
    CNglo   = cell(nDOM,1) ;  % Global connectivities, all nodes
    idom = 1; % Domain 1 
    itype = TYPEslices(idom) ; 
    CNglo{idom} = DATA3D{itype}.CN  ;
    CONNECTbGLO=  cell(nDOM,max(NFACES_REF)) ; % Global connectivities, boundary nodes
    CONNECTbGLO(idom,1:length(CONNECTb_faces{itype})) =  CONNECTb_faces{itype};    
    % We create also an array for storing the nodes of faces F1 and F2 of
    % each domain.
    NODES_faces12 = cell(nDOM,2) ;
    NODES_faces12{idom,1} = f1NOD{itype} ;
    NODES_faces12{idom,2} = f2NOD{itype} ;
    % ---------------------------
    Materials =cell(nDOM,1) ;
    Materials{idom} =DATA3D{itype}.MaterialType;
    MATERIALNEW.PLY = MATERIAL(itype).PLY ; 
    % --------------------------
    ListElementsDom = cell(nDOM,1) ;
    ListElementsDom{idom} = 1:size(DATA3D{itype}.CN,1) ;  % List of elements domain 1
    ListNodesDom = cell(nDOM,1) ;
    ListNodesDom{idom} = 1:size(DATA3D{itype}.COOR,1) ;  % List of nodes domain 1
    
    %% Loop over domains     
    nnodesACUM = size(DATA3D{itype}.COOR,1) ; 
    nelemACUM = size(DATA3D{itype}.CN,1) ; 
    nmatACUM = nmat(itype) ; 
    MaterialTypeNoSlices = DATA3D{itype}.MaterialType ; 
    
    for e = 2:nDOM
        itype = TYPEslices(e) ;
        itype_prev = TYPEslices(e-1) ;
        disp(['Cell =',num2str(e),'; TYPE = ',num2str(itype)])    ;    
        rotGLO{e}  = DATA.ANGiniGLO{e} ;  
        %% List of nodes
        % ----------------------------------------------------
        nnodes = size(DATA3D{itype}.COOR,1) ; 
        nelem = size(DATA3D{itype}.CN,1) ; 
     
        ListNodesDom{e} = (1:nnodes)+nnodesACUM ; % List of nodes domain ( for ensuring same numbering
        % of nodes for all domains)
        % Connectivity matrix domain. All elements
        %----------------------------------------------
        CNnew = DATA3D{itype}.CN + nnodesACUM; % We sum up the number of nodes of each domain
        ListElementsDom{e} = (1:nelem)+nelemACUM;
        % ---------------------------------------------------------------------------------
        % Boundary elements (for each labeled face)
        % ---------------------------------------------------------------------------------
        %
        %     %    CONNECTbGLO(idom,1:length(CONNECTb_faces{itype})) =  CONNECTb_faces{itype};
        %         CNbLOC = CONNECTbGLO(e-1,:) ;
        %         for iface =1:length(CONNECTb_faces)
        %             CNbLOC{iface} = CONNECTbGLO{e-1,iface} + nnodes;
        %         end
        CNbLOC = cell(size(CONNECTb_faces{itype})) ;
        for iface= 1:length(CONNECTb_faces{itype})
            CNbLOC{iface} = CONNECTb_faces{itype}{iface} + nnodesACUM ;
        end
        
        % 
        f1NEW = f1NOD{itype} + nnodesACUM;
        f2NEW = f2NOD{itype} + nnodesACUM;
        f2OLD = NODES_faces12{e-1,2} ;
      
        %         % -------------------------------------------------------------
%       NODES_faces12{e,1} = NODES_faces12{e-1,2} ;   % Face 1,
%       NODES_faces12{e,2} = NODES_faces12{e-1,2} + nnodesACUM ;   % Face 2
        NODES_faces12{e,1} = f2OLD ;   % Face 1, --- Continuity condition 
        NODES_faces12{e,2} = f2NEW ; 
%         % -----------------------------------------------------------------------------------
        % ----------        
        % Renumbering: we have to replace nodes (f1NEW) by f2NOD (this means that f1NEW dissapears from the CN matrices)
        % in CNnew and CNbLOC
        ListNodesDom{e}(f1NOD{itype}) =  ListNodesDom{e-1}(f2NOD{itype_prev}) ;        
        CNnewREN  = CNnew ;
        CNbLOCnew = CNbLOC ;
        for ifacen = 1:length(f1NOD{itype})   % This may be improved (make it more efficient)...
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
        CONNECTbGLO(e,1:length(CNbLOCnew)) = CNbLOCnew ; %
        Materials{e} = DATA3D{itype}.MaterialType + nmatACUM ;
        MaterialTypeNoSlices = [MaterialTypeNoSlices; DATA3D{itype}.MaterialType] ; 
        for imat = 1:length(MATERIAL(itype).PLY)
            matREF = nmatACUM;
            MATERIALNEW.PLY(matREF+imat) =  MATERIAL(itype).PLY(imat) ;
        end      
          nnodesACUM = nnodesACUM + nnodes; 
        nelemACUM = nelemACUM + nelem ; 
        nmatACUM = nmatACUM + nmat(itype) ; 
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
DATA.MaterialTypeNoSlices = MaterialTypeNoSlices; 
if DATA.PREPRINT_MESH == 1
    d = [] ; strainGLOgid = [] ;  stressGLOgid = []; React = [] ; posgp = [] ;
    NameFileMesh = ['prueba.msh'] ;
    Fnodes = [] ;    
    NAME_BASE_GIDfiles = DATA.INPUTDATAfile ; %fileparts() ;     
    GidPostProcess(COOR,CN,TypeElement,d,strainGLOgid, stressGLOgid,  ...
        React,NAME_BASE_GIDfiles,posgp,NameFileMesh,Materials,DATA,Fnodes);     
    error('Disable this error command to go on !!!!')    
end




end

