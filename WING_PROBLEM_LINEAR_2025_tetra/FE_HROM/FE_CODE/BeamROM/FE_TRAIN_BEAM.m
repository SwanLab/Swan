function FE_TRAIN_BEAM(NameFileMeshLOC,NameFileMeshLOC_coarse,NAME_DATA_REF,NAME_ROOT_FE_BEAM,EXECUTABLE_FOLDER,...
    DATAin)

if nargin == 0
    EXECUTABLE_FOLDER = '/home/joaquin/Desktop/CURRENT_TASKS/COMPOSITE_MATERIALS_DOCENCIA/APUNTES_LATEX/DOCUMENTOS_anexos/MATLAB/ELASTOSTATIC_GEN/' ;
    NameFileMeshLOC =   'MESH_slices/TUBE_COMP/FINE_18k' ;  'MESH_slices/TUBE_COMP/FINE_18k' ;
    NameFileMeshLOC_coarse ='MESH_slices/TUBE_COMP/COARSE_1mat_630' ;  'MESH_slices/TUBE_COMP/COARSE_1000' ;
    NAME_ROOT_FE_BEAM = 'TEST_' ;
    NAME_DATA_REF = 'DATAREF' ;
    
elseif nargin == 3
    NAME_ROOT_FE_BEAM = 'TEST_' ;
    EXECUTABLE_FOLDER = '/home/joaquin/Desktop/CURRENT_TASKS/COMPOSITE_MATERIALS_DOCENCIA/APUNTES_LATEX/DOCUMENTOS_anexos/MATLAB/ELASTOSTATIC_GEN/' ;
    
end

% Folder in which the main program is located.
FOLDER = cd ;
TYPE_BOUNDARY_CONDITIONS =  'PERIODIC_BEAMS_MASS_MATRIX';

%DATA.INCLUDE_SELFWEIGTH = 1;
%%%%%%%%%%%%%%%%%
nprojects = 6 ;
IMPOSED_DISP.RIGHT_END = zeros(nprojects,6) ;
IMPOSED_DISP.LEFT_END = zeros(nprojects,6) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ndom = DATAin.NSLICES;

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
