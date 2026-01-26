function  FE_TRAIN_RVE2D(NameFileMeshLOC,NameFileMeshLOC_coarse,...
    NAME_DATA_REF,NAME_ROOT_FE_BEAM,EXECUTABLE_FOLDER,NUMBER_OF_RVES,PRESCRIBED_VALUE,DATARUN)

% Copy of FE_TRAIN_BEAM_COMBINED.m.
% JAHO, 21-July-2018
% --------------------------------------------------------------------------


if nargin == 0
    load('tmp2.mat')
end

% Folder in which the main program is located.
FOLDER = cd ;

% Defining set of FE training tests
% ---------------------------------

[NAMEPROJECTS,IMPOSED_DISP,NUMBER_OF_DOMAINS,NAMES_LOC_PROJ,IMPOSED_FORCE] =...
    NameProjectsTraining2Drve(PRESCRIBED_VALUE,NUMBER_OF_RVES,NAME_ROOT_FE_BEAM,...
    DATARUN) ;







BoundaryConditions = cell(size(NAMEPROJECTS)) ;
%BoundaryConditions(1:4) = {'PERIODIC_BEAMS_MASS_MATRIX'};
switch DATARUN.BoundaryConditionType
    case {'PRESCRIBED_DISP_PLATES'}
        BoundaryConditions(:) = {'PRESCRIBED_DISP_PLATES_RVE'} ;
    case {'PRESCRIBED_DISP_PLATES_FREE'}
        BoundaryConditions(:) = {'PRESCRIBED_DISP_PLATES_RVE_FREE'} ;
    case {'PRESCRIBED_DISPLACEMENTS','PRESCRIBED_DISPLACEMENTS_AND_FORCES','PRESCRIBED_DISPLACEMENTS_BEAMLIKE'}
        BoundaryConditions(:) = {'FIXED_FACES_RVES'} ;
    case 'PRESCRIBED_DISPLACEMENTS_LINEAR'
        BoundaryConditions(:) = {'FIXED_FACES_LINEAR_LATERAL_RVES'} ;
    case 'MIXED_LINEAR_AND_ZERO_PRESCRIBED'
        nprojFIXEd = length(DATARUN.FACES_WITH_NON_ZERO_DISPLACEMENTS)*6 ;
        BoundaryConditions(1:nprojFIXEd) ={'FIXED_FACES_RVES'}  ;
        BoundaryConditions((nprojFIXEd+1):end) ={'FIXED_FACES_LINEAR_LATERAL_RVES'}  ;
    otherwise
        error('Option not implemented')
end



DATA_TYPESOLVER = 0; % Conjugated gradient if = 1.
DATA_niterCONJG = 100000 ; % Number of iterations  conjugated gradient

FOLDER_SAVE = [FOLDER,'/DATAFE_TRAIN/'] ;
if exist(FOLDER_SAVE)==0;     mkdir(FOLDER_SAVE) ; end
DATA.RECALCULATE_STIFFNESS = 1 ;
DATARUN = DefaultField(DATARUN,'ComputeStiffnessOnlyOnce',1) ;
if ~exist('FE_ELASTOSTATIC','file')
    addpath([EXECUTABLE_FOLDER,'FE_CODE']) ;
end
PROJECT_REFERENCE = [] ;
%dBENDING = [] ;

%DATA_FACES = [] ;
nfaces_max = 8 ;

DATARUN = DefaultField(DATARUN,'PostProcessWithNOSLICES',1) ;
DATARUN = DefaultField(DATARUN,'STORE_STIFFNESS',2) ;

for iprojects = 1:length(NAMEPROJECTS)
    disp('------------------------------')
    disp(['PROJECT = ',num2str(iprojects)])
    disp('------------------------------')
    % NAMELOC = [FOLDER,'/',NAMEPROJECTS{iprojects}] ;
    
    if isempty(NAMEPROJECTS{iprojects})
    else
        %  disp(['DATA.RECALCULATE_STIFFNESS=',num2str(DATA.RECALCULATE_STIFFNESS)]) ;
        DATA.nameWORKSPACE =[FOLDER_SAVE,'/',NAMEPROJECTS{iprojects},'.mat'] ;
        DATA.PostProcessWithNOSLICES = DATARUN.PostProcessWithNOSLICES;
        
        
        
        cd(EXECUTABLE_FOLDER) ;
        run(NAME_DATA_REF) ;
        
        DATA.INPUTDATAfile = [FOLDER,'/',NAMEPROJECTS{iprojects}] ;
        FUNinput.INPUTS.DISPLACEMENT_COND_FILE =BoundaryConditions{iprojects} ;
        FUNinput.INPUTS.PLATELOADS.LOAD_UNIFORM = cell(1,nfaces_max) ;
        FUNinput.INPUTS.PLATELOADS.ISLOCAL = cell(1,nfaces_max) ;
        ndom = prod((NUMBER_OF_DOMAINS{iprojects})) ;
        if  isempty(IMPOSED_FORCE) || isempty(IMPOSED_FORCE{iprojects})
            
            FUNinput.INPUTS.PLATELOADS.LOAD_UNIFORM(:) = {zeros(ndom,3) };
            FUNinput.INPUTS.PLATELOADS.ISLOCAL(:) = {zeros(ndom,1) };
        else
            FUNinput = ImposedForcesPlates(FUNinput,IMPOSED_FORCE{iprojects},ndom,...
                NUMBER_OF_DOMAINS{iprojects})  ;
            
        end
        % FE CODE
        DATA.MakeMeshByRepetition.nDOM = NUMBER_OF_DOMAINS{iprojects} ;
        DATA.angROTATION_FACE = DATARUN.angDOM;
        
        if ~isempty(IMPOSED_DISP)
            for ifacee = 1:length(IMPOSED_DISP)
                DISPLACEMENTS_RIGID_BODY{ifacee}  =  IMPOSED_DISP{ifacee}(iprojects,:) ;
            end
        else
            DISPLACEMENTS_RIGID_BODY  =PRESCRIBED_VALUE;
        end
        
        
        FUNinput.INPUTS.NameFileMeshLOC_coarse = NameFileMeshLOC_coarse ;
        
        FUNinput.INPUTS.NameFileMesh = NameFileMeshLOC ;
        FUNinput.INPUTS.DISPLACEMENTS_RIGID_BODY = DISPLACEMENTS_RIGID_BODY ;
        FUNinput.INPUTS.NAMEPROJECT = NAMES_LOC_PROJ{iprojects} ;
        if isempty(PROJECT_REFERENCE)
            PROJECT_REFERENCE = iprojects ;
        end
        DATA.nameWORKSPACE_Kstiff = [FOLDER_SAVE,'/',NAMEPROJECTS{PROJECT_REFERENCE},'.mat'] ;
        
        DATA.TYPESOLVER  =DATA_TYPESOLVER; % Conjugated gradient if = 1.
        DATA.niterCONJG  = DATA_niterCONJG; % Conjugated gradient if = 1.
        FUNinput.INPUTS.MATERIAL =  MATERIAL   ;
        DATA.TOLERANCE_PERIODIC_NODES_GIVEN = DATARUN.TOLERANCE_PERIODIC_NODES_GIVEN ; 
                DATA.STORE_STIFFNESS = DATARUN.STORE_STIFFNESS ; 

        DATA.CalculateNst = 1;
        
        
        
        AAAA = tic ;
        [DATAOUT,DATAIN_OUT ]= FE_ELASTOSTATIC(FUNinput,DATA) ;
        AAAA = toc(AAAA) ;
        disp(['Total time FE training simulation = ',num2str(AAAA),' s'])
        
        if DATARUN.ComputeStiffnessOnlyOnce==1
            DATA.RECALCULATE_STIFFNESS  = 0;
        end
        
        
        %
        cd(FOLDER) ;
        DATA.RECALCULATE_STIFFNESS =0 ;
        
        
    end
end

% NAME_BINARY_SLICE = [FOLDER_SAVE,filesep,NAME_ROOT_FE_BEAM,'global.mat'] ;
% save(NAME_BINARY_SLICE,'DATA_FACES')

