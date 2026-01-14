function RUN_FE(EXECUTABLE_FOLDER,FOLDER,DATA,NAMEPROJECTS,NAME_INPUT_DATA_MATERIAL,IMPOSED_FORCEbnd,...
    NUMBER_OF_DOMAINS,BoundaryConditions,IMPOSED_DISP,PRESCRIBED_VALUE,NameFileMeshLOC_coarse,NameFileMeshLOC,DATARUN)

% ------------------------------------------------------------------
% isurface= 3; %Internal face
% idirection = 1 ;  % Normal direction, always first component
% islices= [10] ;
% VAL = 100 ;
% LOAD_UNIFORM{iproj}{isurface}(islices,idirection) = VAL ;
% ISLOCAL{iproj}{isurface}(islices) = 1 ;
%INPUTS_LOC.BEAMLOADS = DefaultField(INPUTS_LOC.BEAMLOADS,'ISLOCAL',ISLOCAL_DEF) ;

%IMPOSED_FORCE = [] ;

DATA_TYPESOLVER = 0; % Conjugated gradient if = 1.
DATA_niterCONJG = 100000 ; % Number of iterations  conjugated gradient
%DATA.USE_PRECOND = 'CHOLESKY';


FOLDER_SAVE = [FOLDER,'/DATAFE_TRAIN/'] ;
if exist(FOLDER_SAVE)==0;     mkdir(FOLDER_SAVE) ; end

DATA.RECALCULATE_STIFFNESS = 1 ;

if ~exist('FE_ELASTOSTATIC','file')
    addpath([EXECUTABLE_FOLDER,'FE_CODE']) ;
end
if ~exist('FE_NONLINEAR','file')
    addpath([EXECUTABLE_FOLDER,'FE_CODE',filesep,'NONLINEAR_Dyn']) ;
end

PROJECT_REFERENCE = [] ;
nfaces_max = 8 ;


DATAINP = [] ;
DATAINP = DefaultField(DATAINP,'StrainStressWith4Components',1) ;

DATA = DefaultField(DATA,'StrainStressWith4Components',DATAINP.StrainStressWith4Components) ;
ffff= fieldnames(DATAINP) ;
for iii =1:length(ffff)
    locname = ffff{iii} ;
    DATA.(locname) = DATAINP.(locname) ;
end

for iprojects = 1:length(NAMEPROJECTS)
    disp('------------------------------')
    disp(['PROJECT = ',num2str(iprojects)])
    disp('------------------------------')
    % NAMELOC = [FOLDER,'/',NAMEPROJECTS{iprojects}] ;
    
    if isempty(NAMEPROJECTS{iprojects})
    else
        %  disp(['DATA.RECALCULATE_STIFFNESS=',num2str(DATA.RECALCULATE_STIFFNESS)]) ;
        DATA.nameWORKSPACE =[FOLDER_SAVE,'/',NAMEPROJECTS{iprojects},'.mat'] ;
        DATA = DefaultField(DATA,'PostProcessWithNOSLICES',0) ; % =1;
        
        
        
        cd(EXECUTABLE_FOLDER) ;
        run(NAME_INPUT_DATA_MATERIAL) ;
        nfaces_max =  8 ;
        DATA.INPUTDATAfile = [FOLDER,'/',NAMEPROJECTS{iprojects}] ;
        FUNinput.INPUTS.DISPLACEMENT_COND_FILE =BoundaryConditions{iprojects} ;
        % LOADS ON PLATES
        FUNinput.INPUTS.PLATELOADS.LOAD_UNIFORM = cell(1,nfaces_max) ;
        FUNinput.INPUTS.PLATELOADS.ISLOCAL = cell(1,nfaces_max) ;
        ndom = prod((NUMBER_OF_DOMAINS{iprojects})) ;
        if  isempty(IMPOSED_FORCEbnd) || isempty(IMPOSED_FORCEbnd{1})
            
            FUNinput.INPUTS.PLATELOADS.LOAD_UNIFORM(:) = {zeros(ndom,3) };
            FUNinput.INPUTS.PLATELOADS.ISLOCAL(:) = {zeros(ndom,1) };
        else
            FUNinput = ImposedForcesPlatesRVE(FUNinput,IMPOSED_FORCEbnd,ndom,...
                NUMBER_OF_DOMAINS{iprojects})  ;
            
        end
        
        % FE CODE
        DATA.MakeMeshByRepetition.nDOM = NUMBER_OF_DOMAINS{iprojects} ;
        DATARUN = DefaultField(DATARUN,'angDOM',[]) ;
        DATA.angROTATION_FACE = DATARUN.angDOM;
        
        
        for ifacee = 1:length(IMPOSED_DISP)
            DISPLACEMENTS_RIGID_BODY{ifacee}  =  IMPOSED_DISP{ifacee}(iprojects,:) ;
        end
        FUNinput.INPUTS.DISPLACEMENTS_RIGID_BODY = DISPLACEMENTS_RIGID_BODY ;
        
        FUNinput.INPUTS.NameFileMeshLOC_coarse = NameFileMeshLOC_coarse ;
        
        FUNinput.INPUTS.NameFileMesh = NameFileMeshLOC ;
     %   FUNinput.INPUTS.NAMEPROJECT = NAMES_LOC_PROJ{iprojects} ;
        
        % Time evolution of prescribed displacements
        % -------------------------------------------
        PRESCRIBED_VALUE = DefaultField(PRESCRIBED_VALUE,'TIME_FACTOR',[0,1]) ;
        
        
        IS_LINEAR_STATIC  =0 ;
        DATA.FACTOR_TIME_DISPLACEMENTS = PRESCRIBED_VALUE.TIME_FACTOR ;
        PRESCRIBED_VALUE.TIME_INTERVAL = PRESCRIBED_VALUE.TIME_FACTOR ;
        PRESCRIBED_VALUE = DefaultField(PRESCRIBED_VALUE,'TIME_DISCRETIZATION',PRESCRIBED_VALUE.TIME_FACTOR ) ;
        DATA.TIME_DISCRETIZATION = PRESCRIBED_VALUE.TIME_DISCRETIZATION ;
        DATA.FACTOR_TIME_BODY_FORCES = ones(size(PRESCRIBED_VALUE.TIME_FACTOR)) ; % Default, gravity
        DATA.FACTOR_TIME_TRACTION_FORCES =  PRESCRIBED_VALUE.TIME_FACTOR ; %
        
        
        if isempty(PROJECT_REFERENCE)
            PROJECT_REFERENCE = iprojects ;
        end
        DATA.nameWORKSPACE_Kstiff = [FOLDER_SAVE,'/',NAMEPROJECTS{PROJECT_REFERENCE},'.mat'] ;
        
        DATA.TYPESOLVER  =DATA_TYPESOLVER; % Conjugated gradient if = 1.
        DATA.niterCONJG  = DATA_niterCONJG; % Conjugated gradient if = 1.
        FUNinput.INPUTS.MATERIAL =  MATERIAL   ;
        
        DATA.CalculateNst = 1;
        
        
        AAAA = tic ;
        
        FE_NONLINEAR(FUNinput,DATA) ;
        
        AAAA = toc(AAAA) ;
        disp(['Total time FE training simulation = ',num2str(AAAA),' s'])
        
        cd(FOLDER) ;
        DATA.RECALCULATE_STIFFNESS =1 ;  % I disable this to avoid conflicts (27-Sept-2018)
        
    end
end

% NAME_BINARY_SLICE = [FOLDER_SAVE,filesep,NAME_ROOT_FE_BEAM,'global.mat'] ;
% save(NAME_BINARY_SLICE,'DATA_FACES')

