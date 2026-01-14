function  FE_TEST_RVE2D(NameFileMeshLOC,...
    NAME_DATA_REF,NAME_ROOT_FE_BEAM,EXECUTABLE_FOLDER,DATARUN,NUMBER_OF_RVES,DATAINP)

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
NAMEPROJECTS = {NAME_ROOT_FE_BEAM} ;
NUMBER_OF_DOMAINS = {NUMBER_OF_RVES} ;

BoundaryConditions ={DATARUN.BoundaryConditionType};


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

nfaces_max = 8 ;

for iprojects = 1:length(NAMEPROJECTS)
    disp('------------------------------')
    disp(['PROJECT = ',num2str(iprojects)])
    disp('------------------------------')
    % NAMELOC = [FOLDER,'/',NAMEPROJECTS{iprojects}] ;
    
    if isempty(NAMEPROJECTS{iprojects})
    else
        %  disp(['DATA.RECALCULATE_STIFFNESS=',num2str(DATA.RECALCULATE_STIFFNESS)]) ;
        DATA.nameWORKSPACE =[FOLDER_SAVE,'/',NAMEPROJECTS{iprojects},'.mat'] ;
        DATA.PostProcessWithNOSLICES = 1;
        
        
        
        cd(EXECUTABLE_FOLDER) ;
        run(NAME_DATA_REF) ;
        
        DATA.INPUTDATAfile = [FOLDER,'/',NAMEPROJECTS{iprojects}] ;
        FUNinput.INPUTS.DISPLACEMENT_COND_FILE =BoundaryConditions{iprojects} ;
     %   FUNinput.INPUTS.BEAMLOADS.LOAD_UNIFORM = cell(1,nfaces_max) ;
 %       ndom = prod((NUMBER_OF_DOMAINS{iprojects})) ;
 
        FUNinput.INPUTS.PLATELOADS.LOAD_UNIFORM = DATAINP.LOAD_UNIFORM;
        % FE CODE
        DATA.MakeMeshByRepetition.nDOM = NUMBER_OF_DOMAINS{iprojects} ;
        
        %         if ~isempty(IMPOSED_DISP)
        %             for ifacee = 1:length(IMPOSED_DISP)
        %                 DISPLACEMENTS_RIGID_BODY{ifacee}  =  IMPOSED_DISP{ifacee}(iprojects,:) ;
        %             end
        %         else
        DISPLACEMENTS_RIGID_BODY  =DATAINP.IMPOSED_DISP;
        %   end
        
        
         
        FUNinput.INPUTS.NameFileMesh = NameFileMeshLOC ;
        FUNinput.INPUTS.DISPLACEMENTS_RIGID_BODY = DISPLACEMENTS_RIGID_BODY ;
        if isempty(PROJECT_REFERENCE)
            PROJECT_REFERENCE = iprojects ;
        end
        DATA.nameWORKSPACE_Kstiff = [FOLDER_SAVE,'/',NAMEPROJECTS{PROJECT_REFERENCE},'.mat'] ;
        
        DATA.TYPESOLVER  =DATA_TYPESOLVER; % Conjugated gradient if = 1.
        DATA.niterCONJG  = DATA_niterCONJG; % Conjugated gradient if = 1.
        FUNinput.INPUTS.MATERIAL =  MATERIAL   ;
        
        DATA.CalculateNst = 1;
        
        
        
        AAAA = tic ;
        [DATAOUT,DATAIN_OUT ]= FE_ELASTOSTATIC(FUNinput,DATA) ;
        AAAA = toc(AAAA) ;
        disp(['Total time FE training simulation = ',num2str(AAAA),' s'])
        
        
        
        
        %
        cd(FOLDER) ;
        DATA.RECALCULATE_STIFFNESS =0 ;
        
        
    end
end

% NAME_BINARY_SLICE = [FOLDER_SAVE,filesep,NAME_ROOT_FE_BEAM,'global.mat'] ;
% save(NAME_BINARY_SLICE,'DATA_FACES')

