function [NAME_PROJ_LOC ]= FE_TRAIN_BEAM_ALLMIXED(NameFileMeshLOC,NameFileMeshLOC_coarse,...
    NAME_DATA_REF,NAME_ROOT_FE_BEAM,EXECUTABLE_FOLDER,NUMBER_OF_SLICES,PRESCRIBED_VALUE)

% Copy of FE_TRAIN_BEAM.m. Prescribed forces rather than displacements
% JAHO, 10-June-2018
% --------------------------------------------------------------------------


if nargin == 0
    load('tmp.mat')
end

% Folder in which the main program is located.
FOLDER = cd ;
NAME_PROJ_LOC = {'axial','torsion','pbend_y','pbend_z','sbend_y','sbend_z'} ;

%DATA.INCLUDE_SELFWEIGTH = 1;
%%%%%%%%%%%%%%%%%
nprojects = 6 ;
IMPOSED_DISP.RIGHT_END = zeros(nprojects,6) ;
IMPOSED_DISP.LEFT_END = zeros(nprojects,6) ;
VALUE = PRESCRIBED_VALUE.FORCE ; 
VALUE_MOMENT = PRESCRIBED_VALUE.MOMENT ; 

GENERALIZED_FORCES_ENDS_BEAM.RIGHT_END = zeros(nprojects,6) ;  % Generalized forces
GENERALIZED_FORCES_ENDS_BEAM.LEFT_END = zeros(nprojects,6) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ndom = NUMBER_OF_SLICES;

BoundaryConditions = cell(size(NAME_PROJ_LOC)) ;
%BoundaryConditions(1:4) = {'PERIODIC_BEAMS_MASS_MATRIX'};
BoundaryConditions(:) = {'BOUNDARY_COND_MIXED'} ;


% %---------------------------------------
% VALUE = 1e3 ; % Value force
% LENGTH_MOMENT = 2 ;
%%%%%%%%%%%%%%%%%%%%%%%%%%
iproj = 1; % Axial test
NAMEPROJECTS{iproj} = [NAME_ROOT_FE_BEAM,NAME_PROJ_LOC{iproj}];
GENERALIZED_FORCES_ENDS_BEAM.RIGHT_END(iproj,1)= VALUE;
NUMBER_OF_DOMAINS(iproj) = ndom	;
%---------------------------------------
iproj = 2; % Torsion test
NAMEPROJECTS{iproj} = [NAME_ROOT_FE_BEAM,NAME_PROJ_LOC{iproj}];
GENERALIZED_FORCES_ENDS_BEAM.RIGHT_END(iproj,4)= VALUE;
NUMBER_OF_DOMAINS(iproj) = ndom	;
iproj = 3;  % Pure bending test -y
NAMEPROJECTS{iproj} =  [NAME_ROOT_FE_BEAM,NAME_PROJ_LOC{iproj}];;
GENERALIZED_FORCES_ENDS_BEAM.RIGHT_END(iproj,5)= VALUE;
%GENERALIZED_FORCES_ENDS_BEAM.LEFT_END(iproj,5)= -VALUE;
NUMBER_OF_DOMAINS(iproj) = ndom	;
%---------------------------------------
iproj = 4; % Pure bending test -z
NAMEPROJECTS{iproj} =  [NAME_ROOT_FE_BEAM,NAME_PROJ_LOC{iproj}];;
GENERALIZED_FORCES_ENDS_BEAM.RIGHT_END(iproj,6)= VALUE;
%GENERALIZED_FORCES_ENDS_BEAM.LEFT_END(iproj,6)= -VALUE;
NUMBER_OF_DOMAINS(iproj) = ndom	;
%------- Simple bending test y
iproj = 5;
NAMEPROJECTS{iproj} =  [NAME_ROOT_FE_BEAM,NAME_PROJ_LOC{iproj}]; ;
GENERALIZED_FORCES_ENDS_BEAM.RIGHT_END(iproj,3)= VALUE;
GENERALIZED_FORCES_ENDS_BEAM.RIGHT_END(iproj,5)= VALUE_MOMENT;
NUMBER_OF_DOMAINS(iproj) = ndom	;
%----%------- Simple bending test z
iproj = 6;
NAMEPROJECTS{iproj} =  [NAME_ROOT_FE_BEAM,NAME_PROJ_LOC{iproj}]; ;
GENERALIZED_FORCES_ENDS_BEAM.RIGHT_END(iproj,2)= VALUE;
GENERALIZED_FORCES_ENDS_BEAM.RIGHT_END(iproj,6)= VALUE_MOMENT;


NUMBER_OF_DOMAINS(iproj) = ndom	;
%
% %---------------------------------------
% iproj = 7;
% NAMEPROJECTS{iproj} =  [NAME_ROOT_FE_BEAM,NAME_PROJ_LOC{iproj}]; ;
% IMPOSED_DISP.RIGHT_END(iproj,5)= 0.01;
% NUMBER_OF_DOMAINS(iproj) = ndom	;
% %---------------------------------------
% iproj = 8;
% NAMEPROJECTS{iproj} =  [NAME_ROOT_FE_BEAM,NAME_PROJ_LOC{iproj}]; ;
% IMPOSED_DISP.RIGHT_END(iproj,6)= 0.01;
% NUMBER_OF_DOMAINS(iproj) = ndom	;

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
dBENDING = [] ;

DATA_DISPLACEMENT = [] ;

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
        FUNinput.INPUTS.DISPLACEMENT_COND_FILE =BoundaryConditions{iprojects} ;
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
        GENERALIZED_FORCES_ENDS_BEAM_DATA{1} = GENERALIZED_FORCES_ENDS_BEAM.LEFT_END(iprojects,:) ;
        GENERALIZED_FORCES_ENDS_BEAM_DATA{2} = GENERALIZED_FORCES_ENDS_BEAM.RIGHT_END(iprojects,:) ;
        
        FUNinput.INPUTS.GENERALIZED_FORCES_ENDS_BEAM = GENERALIZED_FORCES_ENDS_BEAM_DATA ;
        
        
        AAAA = tic ; 
        [DATAOUT,DATAIN_OUT ]= FE_ELASTOSTATIC(FUNinput,DATA) ;
        AAAA = toc(AAAA) ; 
        disp(['Total time FE training simulation = ',num2str(AAAA),' s'])
        
        % We save the displacement vector arising from the bending tests (iproj= 1, and iproj =2 )
        if iprojects == 3 || iprojects == 4
            dBENDING = [dBENDING, DATAOUT.d] ;
            if iprojects ==4
                % Saving the information in a binary file
                DATA.NameWS_bending_displacements = [FOLDER_SAVE,'/','BendingDisplacements','.mat'] ;
                save(DATA.NameWS_bending_displacements,'dBENDING')
            end
        end
        
        %%%%%%%%
        DATA_DISPLACEMENT_faceB.(NAME_PROJ_LOC{iprojects}).DISP =  DATAOUT.d(DATAIN_OUT.DOFB) ; % Displacement end B
        DATA_DISPLACEMENT_faceB.(NAME_PROJ_LOC{iprojects}).INPUT_FORCE =  GENERALIZED_FORCES_ENDS_BEAM_DATA{2}; % Force applied on face B
        
        %%%%%%%
        
        %
        cd(FOLDER) ;
        DATA.RECALCULATE_STIFFNESS =1 ;
        
        
    end
end

NAME_BINARY_SLICE = [FOLDER_SAVE,filesep,NAME_ROOT_FE_BEAM,'global.mat'] ;
save(NAME_BINARY_SLICE,'DATA_DISPLACEMENT_faceB')

