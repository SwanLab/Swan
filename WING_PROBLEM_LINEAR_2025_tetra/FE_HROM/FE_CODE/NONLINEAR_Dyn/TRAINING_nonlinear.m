function TRAINING_nonlinear(EXECUTABLE_FOLDER,FOLDER,DATA,NAMEPROJECTS,NAME_INPUT_DATA_MATERIAL,LOAD_UNIFORM,ISLOCAL,...
    NUMBER_OF_DOMAINS,BoundaryConditions,IMPOSED_DISP,PRESCRIBED_VALUE,NameFileMeshLOC_coarse,NameFileMeshLOC)






% isurface= 3; %Internal face 
% idirection = 1 ;  % Normal direction, always first component
% islices= [10] ;
% VAL = 100 ;
% LOAD_UNIFORM{iproj}{isurface}(islices,idirection) = VAL ;
% ISLOCAL{iproj}{isurface}(islices) = 1 ;
%INPUTS_LOC.BEAMLOADS = DefaultField(INPUTS_LOC.BEAMLOADS,'ISLOCAL',ISLOCAL_DEF) ;



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
dBENDING = [] ;

DATA_FACES = [] ;
DATAINP = [] ; 
DATAINP = DefaultField(DATAINP,'StrainStressWith4Components',1) ; 
DATA = DefaultField(DATA,'DATAINP',[]) ; 
DATA.DATAINP = DefaultField(DATA.DATAINP,'OnlyShowPlotReactionForces',0) ; 

DATA.OnlyShowPlotReactionForces = DATA.DATAINP.OnlyShowPlotReactionForces ; 

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
        DATA = DefaultField(DATA,'PostProcessWithNOSLICES',1) ; % =1;
        
        
        
        cd(EXECUTABLE_FOLDER) ;
        run(NAME_INPUT_DATA_MATERIAL) ;
        
        DATA.INPUTDATAfile = [FOLDER,'/',NAMEPROJECTS{iprojects}] ;
        FUNinput.INPUTS.DISPLACEMENT_COND_FILE =BoundaryConditions{iprojects} ;
        if exist('LOAD_UNIFORM','var')
             FUNinput.INPUTS.BEAMLOADS.LOAD_UNIFORM = LOAD_UNIFORM{iprojects} ; 
             FUNinput.INPUTS.BEAMLOADS.ISLOCAL = ISLOCAL{iprojects} ; 
        else
        FUNinput.INPUTS.BEAMLOADS.LOAD_UNIFORM = cell(1,5) ;
        FUNinput.INPUTS.BEAMLOADS.LOAD_UNIFORM(:) = {zeros(NUMBER_OF_DOMAINS(iprojects),2) };
        end
        
        
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
        
        % Time evolution of prescribed displacements 
        % -------------------------------------------
        PRESCRIBED_VALUE = DefaultField(PRESCRIBED_VALUE,'TIME_FACTOR',[]) ; 
        
        if isempty(PRESCRIBED_VALUE.TIME_FACTOR)
            % Linear, static problem problem 
            IS_LINEAR_STATIC  =1 ; 
        else
             IS_LINEAR_STATIC  =0 ; 
             DATA.FACTOR_TIME_DISPLACEMENTS = PRESCRIBED_VALUE.TIME_FACTOR ; 
             PRESCRIBED_VALUE.TIME_INTERVAL = PRESCRIBED_VALUE.TIME_DISCRETIZATION ; 
             PRESCRIBED_VALUE = DefaultField(PRESCRIBED_VALUE,'TIME_DISCRETIZATION',PRESCRIBED_VALUE.TIME_FACTOR ) ; 
             DATA.TIME_DISCRETIZATION = PRESCRIBED_VALUE.TIME_DISCRETIZATION ;              
             DATA.FACTOR_TIME_BODY_FORCES = ones(size(PRESCRIBED_VALUE.TIME_FACTOR)) ; % Default, gravity
             DATA.FACTOR_TIME_TRACTION_FORCES =  PRESCRIBED_VALUE.TIME_FACTOR ; % 
             
             
        end
      
        
        if isempty(PROJECT_REFERENCE)
            PROJECT_REFERENCE = iprojects ;
        end
        DATA.nameWORKSPACE_Kstiff = [FOLDER_SAVE,'/',NAMEPROJECTS{PROJECT_REFERENCE},'.mat'] ;
        
        DATA.TYPESOLVER  =DATA_TYPESOLVER; % Conjugated gradient if = 1.
        DATA.niterCONJG  = DATA_niterCONJG; % Conjugated gradient if = 1.
        FUNinput.INPUTS.MATERIAL =  MATERIAL   ;
        
        DATA.CalculateNst = 1;
     %   GENERALIZED_FORCES_ENDS_BEAM_DATA{1} = GENERALIZED_FORCES_ENDS_BEAM.LEFT_END(iprojects,:) ;
     %   GENERALIZED_FORCES_ENDS_BEAM_DATA{2} = GENERALIZED_FORCES_ENDS_BEAM.RIGHT_END(iprojects,:) ;
        
    %    FUNinput.INPUTS.GENERALIZED_FORCES_ENDS_BEAM = GENERALIZED_FORCES_ENDS_BEAM_DATA ;
        
        
        AAAA = tic ; 
        if IS_LINEAR_STATIC == 1
         FE_ELASTOSTATIC(FUNinput,DATA) ;
        else
            DATA = DefaultField(DATA,'OnlyShowPlotReactionForces',0) ; 
            
            if DATA.OnlyShowPlotReactionForces ==0
                [DATAOUT, DATA,OPERfe,NODES_SNAP]=   FE_NONLINEAR(FUNinput,DATA) ;
                
                ReactionsPlot(DATA,OPERfe,NODES_SNAP) ;
                
            else
                
                PlotStoredReactions(DATA.nameWORKSPACE) ;
                
            end
            
            
        end
        AAAA = toc(AAAA) ; 
        
        
        
        disp(['Total time FE training simulation = ',num2str(AAAA),' s'])
        
        % We save the displacement vector arising from the bending tests (iproj= 1, and iproj =2 )
%         if iprojects == 3 || iprojects == 4
%             dBENDING = [dBENDING, DATAOUT.d] ;
%             if iprojects ==4
%                 % Saving the information in a binary file
%                 DATA.NameWS_bending_displacements = [FOLDER_SAVE,'/','BendingDisplacements','.mat'] ;
%                 save(DATA.NameWS_bending_displacements,'dBENDING')
%             end
%         end
        
        %%%%%%%%
    %    DATA_FACES.(NAME_PROJ_LOC{iprojects}).DISP =  DATAOUT.d(DATAIN_OUT.DOFB) ; % Displacement end B
     %   DATA_FACES.(NAME_PROJ_LOC{iprojects}).REACT =  DATAOUT.React(DATAIN_OUT.DOFB) ; % Displacement end B
        
        
    %    DATA_DISPLACEMENT_faceB.(NAME_PROJ_LOC{iprojects}).INPUT_FORCE =  GENERALIZED_FORCES_ENDS_BEAM_DATA{2}; % Force applied on face B
        
        %%%%%%%
        
        %
        cd(FOLDER) ;
        DATA.RECALCULATE_STIFFNESS =1 ;  % I disable this to avoid conflicts (27-Sept-2018)
        
        
    end
end

% NAME_BINARY_SLICE = [FOLDER_SAVE,filesep,NAME_ROOT_FE_BEAM,'global.mat'] ;
% save(NAME_BINARY_SLICE,'DATA_FACES')

