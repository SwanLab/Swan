function [COOR,CN,TypeElement,TypeElementB, PROPMAT,  DOFr,dR,...
    Tnod,CNb,fNOD,Fpnt,NameFileMesh,typePROBLEM,celasgloINV,CONNECTb,DOFm,Gb,...
    DATA,MaterialType] = INPUTS_GENERAL_NONLINEAR(INPUTS_LOC,DATA) ;

% Inputs   laminate
%dbstop('7')
if nargin == 0
    load('tmp1.mat')
end


DOFm = [] ; Gb = [] ; DOFA = [] ; DOFB = [] ;
Fpnt_Dirichlet = [] ; Tnod = [] ; CNb = [] ;
DATAcoarsefine = [] ;

% ---------------
% 1.  Finite element mesh:  COORDINATES AND CONNECTIVITIES for both the volume domain and the boundary domain
% OUTPUT: COOR,CN,TypeElement,CONNECTb,TypeElementB

NameFileMesh = INPUTS_LOC.NameFileMesh ;
READ_MATERIAL_COLUMN = 1;

[dummy1 NameFileMeshHERE dummy2]= fileparts(NameFileMesh) ;


nameWORKSPACE = ['DATAWS/',DATA.INPUTDATAfile,'_WS.mat'] ; % To retrieve already read data

DATA = DefaultField(DATA,'nameWORKSPACE',nameWORKSPACE);


DATA= DefaultField(DATA,'FactorDivideCoordinates',1) ;
DATA= DefaultField(DATA,'nameWORKSPACE_Kstiff',[]) ;
DATA.MATERIAL_ORIGINAL = INPUTS_LOC.MATERIAL ;
DATA.DOMAINVAR = [] ;
DOMAINVAR= [] ;

if DATA.RECALCULATE_STIFFNESS == 1
    disp('Reading GID mesh ...')
    DATA = DefaultField(DATA,'MakeMeshByRepetition',[]) ;
    DATA.MakeMeshByRepetition = DefaultField(DATA.MakeMeshByRepetition,'nDOM',[]) ;
    DATA.MakeMeshByRepetition = DefaultField(DATA.MakeMeshByRepetition,'DIRECTION','z') ;
    if isempty(DATA.MakeMeshByRepetition.nDOM)
        [COOR,CN,TypeElement,CONNECTb,TypeElementB,MaterialType]=...
            ReadMeshFile(NameFileMesh,'READ_MATERIAL_COLUMN',READ_MATERIAL_COLUMN)  ;
    else
        DATA = DefaultField(DATA,'FACES_LABELLED_WITH_GID',0) ;  %
        if  DATA.FACES_LABELLED_WITH_GID == 0
            % OLD VERSION OF THE PROGRAM (Before 6-April-2018), Set
            % DATA.FACES_LABELLED_WITH_GID = 0 in the input data file to
            % enable this option
            [COOR,CN,TypeElement,CONNECTb,TypeElementB,MaterialType,INPUTS_LOC.MATERIAL]=...
                MeshGenerationRepeat3D(NameFileMesh,DATA,INPUTS_LOC.MATERIAL);
        else
            % It requires defining faces of the slice via GID problemtype
            if  length(DATA.MakeMeshByRepetition.nDOM) == 1
                % Tiling copies in 1 direction
                [COOR,CN,TypeElement,CONNECTb,TypeElementB,MaterialType,INPUTS_LOC.MATERIAL,DOMAINVAR]=...
                    MeshGenerationRepeat3D_faces(NameFileMesh,DATA,INPUTS_LOC.MATERIAL);
            elseif length(DATA.MakeMeshByRepetition.nDOM) == 2
                % Tiling copies in 2 direction
                [COOR,CN,TypeElement,CONNECTb,TypeElementB,MaterialType,INPUTS_LOC.MATERIAL,DOMAINVAR]=...
                    MeshGenerationRepeat3D_2repeat(NameFileMesh,DATA,INPUTS_LOC.MATERIAL);
            else
                error('Option not implemented')
            end
        end
    end
    COOR = COOR/DATA.FactorDivideCoordinates ;
    [FFOLDER FFILE] =   fileparts(DATA.nameWORKSPACE) ;
    if ~exist(FFOLDER,'dir')
        mkdir(FFOLDER)
    end
    % Renumbering
    %DATA.Renumber = 0 ;
    % if DATA.Renumber == 1
    DATA.RENUMBERED_OUTSIDE = 1 ; % To avoid performing this operation later, when computing stiffness matrix
    [~,IndicesRenumberingElements]  = sort(CN(:,1)) ;
    CN = CN(IndicesRenumberingElements,:) ;
    if ~isempty(MaterialType)
        MaterialType = MaterialType(IndicesRenumberingElements) ;
    end
    %    else
    %        IndicesRenumberingElements = 1:length(MaterialType) ;
    %        IndicesRenumberingElements = IndicesRenumberingElements(:) ;
    %    end
    
    save(DATA.nameWORKSPACE,'COOR','CN','TypeElement','CONNECTb','TypeElementB',...
        'MaterialType','NameFileMesh','DOMAINVAR','IndicesRenumberingElements')
    
    disp('Done')
else
    disp('Retrieving mesh data... DONE')
    NAMEWS = DATA.nameWORKSPACE_Kstiff ;
    if isempty(NAMEWS)
        NAMEWS = DATA.nameWORKSPACE ;
    end
    load(NAMEWS,'COOR','CN','TypeElement','CONNECTb','TypeElementB','MaterialType','DOMAINVAR') ;
    % To avoid problems ....
    save(DATA.nameWORKSPACE,'COOR','CN','TypeElement','CONNECTb','TypeElementB',...
        'MaterialType','NameFileMesh','DOMAINVAR')
end



% % READING .dat files (faces and lines). Generated with GID's problemtype PROBLEMTYPES_GID/SLICE_JOINTS.gid
% nameDAT1 = [NameFileMesh,'.dat'] ;
% if exist()
% [NODES_FACES,NODES_LINES] = NodesFacesLinesGID(NameFileMesh) ;



nnode = size(COOR,1) ;% Number of nodes
ndim = size(COOR,2); % Number of spatial dimensions (ndim=2 for 2D problems)
nelem = size(CN,1) ; % Number of elements

%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. MATERIAL PROPERTIES: output celasglo   %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATERIAL.FIBER.Eyoung = 80000 ; %MN/m2
% MATERIAL.FIBER.POISSON = 0.25 ;
% MATERIAL.FIBER.INDEX = 1;
% %%%
% MATERIAL.MATRIX.Eyoung = 20000    ; 200e7 ; %N/m2
% MATERIAL.MATRIX.POISSON = 0.25 ;
% MATERIAL.MATRIX.INDEX = 2;
%dbstop('125')
MATERIAL = INPUTS_LOC.MATERIAL  ;
celasgloINV = [] ;
[PROPMAT, typePROBLEM, densGLO,DATA,Cglo]  = ...
    AssignMatProp_J2(ndim,MATERIAL,nelem,MaterialType,DATA,TypeElement,size(CN,2)) ;
save(DATA.nameWORKSPACE,'Cglo','-append')

%[celasglo, celasgloINV, typePROBLEM, densGLO,DATA]  = AssignMatPropLAM(ndim,MATERIAL,nelem,MaterialType,DATA) ;

DATA.MATERIAL = MATERIAL ;
DATA.densGLO = densGLO ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Dirichlet (essential) boundary conditions, OUTPUT: dR and rdof
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%dbstop('85')
INPUTS_LOC = DefaultField(INPUTS_LOC,'NameFileMeshLOC_coarse',INPUTS_LOC.NameFileMesh) ;
NameFileMeshLOC_coarse = INPUTS_LOC.NameFileMeshLOC_coarse ;
if exist(DATA.nameWORKSPACE,'file')
    APPEND = '-append' ;
else
    APPEND = '' ;
    
end
save(DATA.nameWORKSPACE,'NameFileMeshLOC_coarse',APPEND) ;

AREA  = [];
NODES_ENTITIES = [] ; BasisUrb_ENTITIES = [] ; 
% CALLING FILE WHICH IS USED TO DEFINE DIRICHLET BOUNDARY CONDITIONS
disp('-------------------------------------------------------') ; 
disp(['DIRICH. BND COND using FILE = ',INPUTS_LOC.DISPLACEMENT_COND_FILE]) ; 
disp('-------------------------------------------------------') ;  
eval(INPUTS_LOC.DISPLACEMENT_COND_FILE);
DATA.AREA = AREA ;  % Area cross-section 
DATA.NODES_ENTITIES  = NODES_ENTITIES ;  % List of nodes of each entity ---boundary
DATA.BasisUrb_ENTITIES = BasisUrb_ENTITIES; % Rigid body modes of each entity. Use for computing reactions resultant

DATA.DOMAINVAR = DOMAINVAR;
DATA.DOFA = DOFA ; % For beams, DOFs initial face
DATA.DOFB = DOFB ; %                  end face
DATA.DATAcoarsefine = DATAcoarsefine;
%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. Neumann (natural) boundary conditions : OUTPUT: Tnod, CNb, Fnod
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISTRIBUTED LOADS and POINT LOADS
% ------------------------
%dbstop('115')
%DATA.NODES_ENTITIES = [] ;
eval(INPUTS_LOC.NEUMANN_BOUNDARY_CONDITIONS);

%%%%%%%%%%%%%%%%%%%%%%5


% POINT LOADS
% -----------
if ~isempty(Fpnt_Dirichlet)  % Nodal forces calculated when computing Dirichlet BCs
    Fpnt =  Fpnt_Dirichlet ;
end

% Body forces....
fNOD = zeros(nnode*ndim,1) ;

%fNOD(933*3) = 0.1 ; % Mpa


%%% POST-PROCESS  and other options
DATA.PLOT.REACTIONS = 1 ;
DATA.CALCULATE_averageSTRESS = 0 ; % Generalized stresses
DATA.VECTcode  = 1 ;   % Vectorize code
DATA.BCBformulation = 1 ;  % Formulation B^T C B
DATA.NOVOIDS = 1 ;
DATA.plotFLUCT  =0;
%DATA.dispMACRO = dispMACRO;
%DATA.NODESpl = NODESpl ;