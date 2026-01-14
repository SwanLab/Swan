function FE_TRAIN_JOINTS(NameFileMeshLOC,NameFileMeshLOC_coarse,NAME_DATA_REF,...
    NAME_ROOT_FE_BEAM,EXECUTABLE_FOLDER,MESH3D,MESH1D,ijoint)

if nargin == 0
   load('tmp.mat')
end

% Folder in which the main program is located.
FOLDER = cd ;
TYPE_BOUNDARY_CONDITIONS =  'PRESCRIBED_FLUCTUATIONS_JOINTS_2faces';  % Fluctuations are given 

%%% Determining fluctuation matrix  (preparing data for FE solver)
% ----------------------------------------------------------------
% How many faces have the joint under consideration  ? 
% First we have to determine which 1D elements are joint 
AREJOINTS = find(MESH1D.TYPE.ISBEAM ==0); 
% Which elements have index ijoint (either joint or beam) ? 
ISINDEX = find(MESH1D.TYPE.INDEX3D ==ijoint); 
ELEM_JOINTS = intersect(AREJOINTS,ISINDEX) ; 
if isempty(ELEM_JOINTS)
    error('This joint does not exist')
end
% We take the first joint 
INDEX = ELEM_JOINTS(1) ; 
% Finally, the number of faces can be consulted in 
N_FACES = MESH1D.NNODES_elem(INDEX) ; 
if N_FACES >2 
    error('This routine only works for 2-faces joint')
end
% And which are the slices connected to the joint ? 
SLICES_CONNECTED = MESH1D.CN(INDEX,:) ; 
% Which 3D indices ? 
INDEX_SLICES_LOC = MESH1D.ELEMENTS_MESH1D(SLICES_CONNECTED) ; 
INDEX_MESH_SLICES = MESH1D.TYPE.INDEX3D(INDEX_SLICES_LOC) ; 

% Now we extract the fluctuation information of each slice. 
% Since we are dealing with straight structures, there is no need to worry
% about rotation matrices  
 FLUCTUATION_MATRIX = cell(1,2) ; 
  RIGID_BODY_MATRIX = cell(1,2) ;
  MASS_MATRIX_GEOMETRIC = cell(1,2) ;
for iii = 1:length(INDEX_MESH_SLICES)
    IND_SLICE = INDEX_MESH_SLICES(iii) ; 
    % Name binary file containing fluctuations modes
    NAMEWSloc = ['MODES',filesep,'MODES_',MESH3D.SLICES(IND_SLICE).NameWSmodes,'.mat'] ; 
    load(NAMEWSloc,'DATAROM','DATA_REFMESH'); 
     FLUCTUATION_MATRIX{iii} = DATAROM.BasisINTfluct ; 
     RIGID_BODY_MATRIX{iii} = DATAROM.BasisIntRB; 
     MASS_MATRIX_GEOMETRIC{iii} = DATA_REFMESH.GeometricMassMatrixInterface ; 
end
% Now we save this information in a binary file 
FLUCTUATIONS_BOUNDARY_NameWS = [FOLDER,filesep,'MODES',filesep,'FLUCTUATIONS.mat'] ; 
save(FLUCTUATIONS_BOUNDARY_NameWS,'FLUCTUATION_MATRIX','RIGID_BODY_MATRIX','MASS_MATRIX_GEOMETRIC')  
DATA.FLUCTUATIONS_BOUNDARY_NameWS = FLUCTUATIONS_BOUNDARY_NameWS ; 
% And what about the mesh ???? Rather than reading the mesh from scratch,
% in this case, mesh data should be an input of the problem. 
% We provide this information in yet another binary file 
NAME_WS_MESH = [FOLDER,filesep,'MODES',filesep,'MESH.mat'] ; 

DATA.NAME_WS_MESH_GIVEN = NAME_WS_MESH ; 
DATA3D_given = MESH3D.JOINTS(ijoint) ; 
save(NAME_WS_MESH,'DATA3D_given') ; 


%

 
nprojects = 6 ;
IMPOSED_DISP.RIGHT_END = zeros(nprojects,6) ;
IMPOSED_DISP.LEFT_END = zeros(nprojects,6) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ndom = 1;

NAME_PROJ_LOC = {'axial','torsion','pbend_y','pbend_z','sbend_y','sbend_z'} ;

iproj = 1;
NAMEPROJECTS{iproj} = [NAME_ROOT_FE_BEAM,NAME_PROJ_LOC{iproj}];
IMPOSED_DISP.RIGHT_END(iproj,1)= 0.01;
NUMBER_OF_DOMAINS(iproj) = ndom	;
%---------------------------------------
iproj = 2;
NAMEPROJECTS{iproj} = [NAME_ROOT_FE_BEAM,NAME_PROJ_LOC{iproj}];
IMPOSED_DISP.RIGHT_END(iproj,4)= 0.01;
NUMBER_OF_DOMAINS(iproj) = ndom	;
%---------------------------------------
iproj = 3;
NAMEPROJECTS{iproj} =  [NAME_ROOT_FE_BEAM,NAME_PROJ_LOC{iproj}];;
IMPOSED_DISP.RIGHT_END(iproj,5)= 0.01;
IMPOSED_DISP.LEFT_END(iproj,5)= -0.01;
NUMBER_OF_DOMAINS(iproj) = ndom	;
%---------------------------------------
iproj = 4;
NAMEPROJECTS{iproj} =  [NAME_ROOT_FE_BEAM,NAME_PROJ_LOC{iproj}];;
IMPOSED_DISP.RIGHT_END(iproj,6)= 0.01;
IMPOSED_DISP.LEFT_END(iproj,6)= -0.01;
NUMBER_OF_DOMAINS(iproj) = ndom	;
%---------------------------------------
iproj = 5;
NAMEPROJECTS{iproj} =  [NAME_ROOT_FE_BEAM,NAME_PROJ_LOC{iproj}]; ;
IMPOSED_DISP.RIGHT_END(iproj,5)= 0.01;
NUMBER_OF_DOMAINS(iproj) = ndom	;
%---------------------------------------
iproj = 6;
NAMEPROJECTS{iproj} =  [NAME_ROOT_FE_BEAM,NAME_PROJ_LOC{iproj}]; ;
IMPOSED_DISP.RIGHT_END(iproj,6)= 0.01;
NUMBER_OF_DOMAINS(iproj) = ndom	;

%

% DILATATION_SLICES_x = linspace(1,1,NUMBER_OF_DOMAINS) ; % Dilatation of slices in the x-direction
%DILATATION_SLICES_x = ones(1,NUMBER_OF_DOMAINS);


DATA_TYPESOLVER = 0; % Conjugated gradient if = 1.
DATA_niterCONJG = 100000 ; % Number of iterations  conjugated gradient
%DATA.USE_PRECOND = 'CHOLESKY';
%NAMEPROJECTS = {'TEST_x','TEST_y','TEST_z','TEST_gx','TEST_gy','TEST_gz'} ;
% NAMEPROJECTS = cell(1,6) ;
% NAMEPROJECTS{2} = 'TEST_yp' ;
%
% NAMEPROJECTS{6} = 'TEST_gzp' ;

FOLDER_SAVE = [FOLDER,'/DATAFE_TRAIN/'] ;
if exist(FOLDER_SAVE)==0;     mkdir(FOLDER_SAVE) ; end

DATA.RECALCULATE_STIFFNESS = 1 ;

if ~exist('FE_ELASTOSTATIC','file')
    addpath([EXECUTABLE_FOLDER,'FE_CODE']) ;
end
PROJECT_REFERENCE = [] ;
for iprojects = 1:length(NAMEPROJECTS)
    disp('------------------------------')
    disp(['PROJECT = ',num2str(iprojects)])
    disp('------------------------------')
    % NAMELOC = [FOLDER,'/',NAMEPROJECTS{iprojects}] ;
    
    if isempty(NAMEPROJECTS{iprojects})
    else
        %  disp(['DATA.RECALCULATE_STIFFNESS=',num2str(DATA.RECALCULATE_STIFFNESS)]) ;
        DATA.nameWORKSPACE =[FOLDER_SAVE,'/',NAMEPROJECTS{iprojects},'.mat'] ;
        DATA.PostProcessWithNOSLICES = 0;
        
       
        
        cd(EXECUTABLE_FOLDER) ;
                 run(NAME_DATA_REF) ;

        DATA.INPUTDATAfile = [FOLDER,'/',NAMEPROJECTS{iprojects}] ;
        FUNinput.INPUTS.DISPLACEMENT_COND_FILE =TYPE_BOUNDARY_CONDITIONS ;
        FUNinput.INPUTS.BEAMLOADS.LOAD_UNIFORM = cell(1,5) ;
        FUNinput.INPUTS.BEAMLOADS.LOAD_UNIFORM(:) = {zeros(NUMBER_OF_DOMAINS(iprojects),3) };
        % FE CODE
        DATA.MakeMeshByRepetition.nDOM = NUMBER_OF_DOMAINS(iprojects) ;
        %   DATA.MakeMeshByRepetition.DILATATION_SLICES_x = DILATATION_SLICES_x ;
        %         MATERIAL.PLY(1).E = E ;
        %         MATERIAL.PLY(1).nu = nu ;
        %         MATERIAL.PLY(1).G = G ;
        %         MATERIAL.PLY(1).typePROBLEM = typePROBLEM ;
        DISPLACEMENTS_ENDS_RIGID_BODY{1} =  IMPOSED_DISP.LEFT_END(iprojects,:) ;
        DISPLACEMENTS_ENDS_RIGID_BODY{2} = IMPOSED_DISP.RIGHT_END(iprojects,:) ;
        %         DISPLACEMENTS_ENDS_RIGID_BODY{1}(:) = mat2cell(mat2cell(IMPOSED_DISP.LEFT_END)) ;
        %         DISPLACEMENTS_ENDS_RIGID_BODY{2}(:) =mat2cell(IMPOSED_DISP.LEFT_END);
        %        DISPLACEMENTS_ENDS_RIGID_BODY{2}{IMPOSED_DISP_DIRECT(iprojects)} = IMPOSED_MOVEMENT ;
        
        
        FUNinput.INPUTS.NameFileMeshLOC_coarse = NameFileMeshLOC_coarse ;
        
        FUNinput.INPUTS.NameFileMesh = NameFileMeshLOC ;
        FUNinput.INPUTS.DISPLACEMENTS_ENDS_RIGID_BODY = DISPLACEMENTS_ENDS_RIGID_BODY ;
        
        if isempty(PROJECT_REFERENCE)
            PROJECT_REFERENCE = iprojects ;
        end
        DATA.nameWORKSPACE_Kstiff = [FOLDER_SAVE,'/',NAMEPROJECTS{PROJECT_REFERENCE},'.mat'] ;
        
        DATA.TYPESOLVER  =DATA_TYPESOLVER; % Conjugated gradient if = 1.
        DATA.niterCONJG  = DATA_niterCONJG; % Conjugated gradient if = 1.
        FUNinput.INPUTS.MATERIAL =  MATERIAL   ;
        
        DATA.CalculateNst = 1;
        
        FE_ELASTOSTATIC(FUNinput,DATA) ;
        
        %
        cd(FOLDER) ;
        DATA.RECALCULATE_STIFFNESS = 0 ;
        
        
    end
end
