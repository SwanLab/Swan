function [COOR,CN,TypeElement,TypeElementB, celasglo,  DOFr,dR,...
    Tnod,CNb,fNOD,Fpnt,NameFileMesh,typePROBLEM,celasgloINV,CONNECTb,DOFm,Gb,...
    DATA,MaterialType] = INPUTS_GENERAL(INPUTS_LOC,DATA) ;

% Inputs   laminate
%dbstop('7')
if nargin == 0
    load('tmp3.mat')
end


DOFm = [] ; Gb = [] ; DOFA = [] ; DOFB = [] ;
Fpnt_Dirichlet = [] ;
DATAcoarsefine = [] ;
% SOLVER
% DATA.niterCONJG = 1000 ; % Number of iterations
% DATA.tolCONJG = 1e-6 ; % Tolerance solver
%%%%%%%%%%%%%%%%%
% ---------------
% 1.  Finite element mesh:  COORDINATES AND CONNECTIVITIES for both the volume domain and the boundary domain
% OUTPUT: COOR,CN,TypeElement,CONNECTb,TypeElementB


NameFileMesh = INPUTS_LOC.NameFileMesh ;
DATA = DefaultField(DATA,'READ_MATERIAL_COLUMN',1) ;%  = 1;


[dummy1 NameFileMeshHERE dummy2]= fileparts(NameFileMesh) ;


% DATA = DefaultField(DATA,'INPUTDATAfile',NameFileMeshHERE)  ;
% [PPATH,FFILE] = fileparts(DATA.INPUTDATAfile) ;
% switch PPATH
%     case {'DATA_input',''}

%DATA = DefaultField(DATA ,'INPUTDATAfile',NameFileMesh) ;

nameWORKSPACE = ['DATAWS/',DATA.INPUTDATAfile,'_WS.mat'] ; % To retrieve already read data
%     otherwise
%         FOLDER_SAVE =[PPATH,'/DATAWS/',] ;
%         if exist(FOLDER_SAVE)==0
%             mkdir(FOLDER_SAVE)
%         end
%         DATA.nameWORKSPACE = [FOLDER_SAVE,'/',FFILE,'_WS.mat'] ; % To retrieve already read data
%
% end

DATA = DefaultField(DATA,'nameWORKSPACE',nameWORKSPACE);

%dbstop('31')
% if exist(DATA.nameWORKSPACE)~=2
%     DATA.RECALCULATE_STIFFNESS = 1 ;
% end

DATA= DefaultField(DATA,'FactorDivideCoordinates',1) ;
DATA= DefaultField(DATA,'nameWORKSPACE_Kstiff',[]) ;
DATA.MATERIAL_ORIGINAL = INPUTS_LOC.MATERIAL ;
DATA.DOMAINVAR = [] ;
DOMAINVAR = [] ;

if DATA.RECALCULATE_STIFFNESS == 1
    disp('Reading GID mesh ...')
    
    DATA = DefaultField(DATA,'MakeMeshByRepetition',[]) ;
    DATA.MakeMeshByRepetition = DefaultField(DATA.MakeMeshByRepetition,'nDOM',[]) ;
    DATA.MakeMeshByRepetition = DefaultField(DATA.MakeMeshByRepetition,'DIRECTION','z') ;
    
    
    if isempty(DATA.MakeMeshByRepetition.nDOM)
        [COOR,CN,TypeElement,CONNECTb,TypeElementB,MaterialType]=...
            ReadMeshFile(NameFileMesh,'READ_MATERIAL_COLUMN',DATA.READ_MATERIAL_COLUMN)  ;
        
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
                DATA = DefaultField(DATA,'angROTATION_FACE',[]) ;
                if  length(DATA.angROTATION_FACE) == 1 || length(DATA.angROTATION_FACE) == 0
                    % Version before 8-oct-2020
                    [COOR,CN,TypeElement,CONNECTb,TypeElementB,MaterialType,INPUTS_LOC.MATERIAL,DOMAINVAR]=...
                        MeshGenerationRepeat3D_faces(NameFileMesh,DATA,INPUTS_LOC.MATERIAL);
                else
                    % Version able to deal with slices of varying
                    % curvature (8-Oct-2020)
                    
                    % CONSTANT CURVATURE FOR EACH SLICE
                    DATA = DefaultField(DATA,'COORtransf',[]) ;
                    
                    if isempty(DATA.COORtransf)
                        
                        [COOR,CN,TypeElement,CONNECTb,TypeElementB,MaterialType,INPUTS_LOC.MATERIAL,DOMAINVAR]=...
                            MeshGenerationRepeat3D_facesVCURV(NameFileMesh,DATA,INPUTS_LOC.MATERIAL);
                        
                    else
                        % General method, curve constructed by  cubic interpolation
                        % See NEW_IDEAS.pdf, BEAM with variable
                        % curvature
                        DATA = DefaultField(DATA,'NAME_MESHES_SLICES',[]) ;
                        
                        if isempty(DATA.NAME_MESHES_SLICES)
                            % Just one reference mesh
                            [COOR,CN,TypeElement,CONNECTb,TypeElementB,MaterialType,INPUTS_LOC.MATERIAL,DOMAINVAR]=...
                                MeshGenerationRepeat3D_CUBIC(NameFileMesh,DATA,INPUTS_LOC.MATERIAL);
                        else
                            % Two or more reference meshes (1-Nov-2020, plane from SOFIA to BCN)
                            [COOR,CN,TypeElement,CONNECTb,TypeElementB,MaterialType,INPUTS_LOC.MATERIAL,DOMAINVAR,DATA]=...
                                MeshGenerationRepeat3D_CUBICDIFF(DATA.NAME_MESHES_SLICES,DATA,INPUTS_LOC.MATERIAL);
                            
                        end
                        
                        
                        
                        
                    end
                end
                
                
                
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
    
    save(DATA.nameWORKSPACE,'COOR','CN','TypeElement','CONNECTb','TypeElementB',...
        'MaterialType','NameFileMesh','DOMAINVAR','-v7.3')
    
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



nnode = size(COOR,1) ;% Number of nodes
ndim = size(COOR,2); % Number of spatial dimensions (ndim=2 for 2D problems)
nelem = size(CN,1) ; % Number of elements

%%% RECALCULATE COORDINATES SO THAT THE MIDPLANE --> Z= 0
DATA = DefaultField(DATA,'CENTRATE_MIDPLANE',0) ;

if DATA.CENTRATE_MIDPLANE == 1
    zmin = min(COOR(:,3)) ; zmin = zmin(1) ;
    zmax = max(COOR(:,3)) ; zmax = zmax(1) ;
    zmed = 0.5*(zmin+zmax) ;
    COOR(:,3) = COOR(:,3)-zmed ;
end

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
%DATA.typePROBLEM = typePROBLEM ;
[celasglo, celasgloINV, typePROBLEM, densGLO,DATA]  = AssignMatPropLAM(ndim,MATERIAL,nelem,MaterialType,DATA) ;
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
    APPEND = '-v7.3' ;
    
end
save(DATA.nameWORKSPACE,'NameFileMeshLOC_coarse',APPEND) ;

AREA  = [];
disp('*****************************************************************************')
disp(['Prescribing displacements with file: ',INPUTS_LOC.DISPLACEMENT_COND_FILE])
disp('*****************************************************************************')
BasisUrbA = [] ;
BasisUrbB = [] ;
RotMatrixA = [] ;
RotMatrixB = [] ;
ANG_ROTATION_TOTAL=[];
BasisUrb = [] ;
NODES_ENTITIES = [] ;
eval(INPUTS_LOC.DISPLACEMENT_COND_FILE);   % Typically, FIXED_FACES_RVES...
DATA.AREA = AREA ;  % Area cross-section

DATA.DOMAINVAR = DOMAINVAR;
DATA.DOFA = DOFA ; % For beams, DOFs initial face
DATA.DOFB = DOFB ; %                  end face

DATA.BasisUrbA = BasisUrbA ;
DATA.BasisUrbB = BasisUrbB ;
DATA.RotMatrixA = RotMatrixA ;
DATA.RotMatrixB = RotMatrixB ;
DATA.BasisUrb = BasisUrb;
DATA.DATAcoarsefine = DATAcoarsefine;
DATA.ANG_ROTATION_TOTAL = ANG_ROTATION_TOTAL ;
DATA.NODES_ENTITIES = NODES_ENTITIES ;




%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. Neumann (natural) boundary conditions : OUTPUT: Tnod, CNb, Fnod
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISTRIBUTED LOADS and POINT LOADS
% ------------------------
%dbstop('115')
INPUTS_LOC = DefaultField(INPUTS_LOC,'NEUMANN_BOUNDARY_CONDITIONS','UNIFORM_LOAD_PLATE_STRUCTURES')  ;
disp('*****************************************************************************')
disp(['Prescribing FORCES with file: ',INPUTS_LOC.NEUMANN_BOUNDARY_CONDITIONS])
disp('*****************************************************************************')

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