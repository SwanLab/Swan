function BeamFEtraining(EXECUTABLE_FOLDER,FOLDER,NameFileMeshLOC,NameFileMeshLOC_coarse,NAME_DATA_REF,...
    DATA,NAME_PROJECTS,NAMEFILE)


TYPE_BOUNDARY_CONDITIONS =  'PERIODIC_BEAMS_MASS_MATRIX';
nprojects = 6 ;
IMPOSED_DISP.RIGHT_END = zeros(nprojects,6) ;
IMPOSED_DISP.LEFT_END = zeros(nprojects,6) ;

DATA_FE_WS = cell(1,nprojects) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ndom = 1;
NAMES_PROJECTS_PP = {'pbend_y','pbend_z','axial','tors','sbend_y','sbend_z'} ;
% -------------------------------
%% FIRST THE PURE BENDING TESTS
% -------------------------------
iproj = 1;
NAMEPROJECTS{iproj} = [NAME_PROJECTS,'_',NAMES_PROJECTS_PP{iproj}] ;
IMPOSED_DISP.RIGHT_END(iproj,5)= 0.01;
IMPOSED_DISP.LEFT_END(iproj,5)= -0.01;
NUMBER_OF_DOMAINS(iproj) = ndom	;
%---------------------------------------
iproj = 2;
NAMEPROJECTS{iproj} = [NAME_PROJECTS,'_',NAMES_PROJECTS_PP{iproj}] ;
IMPOSED_DISP.RIGHT_END(iproj,6)= 0.01;
IMPOSED_DISP.LEFT_END(iproj,6)= -0.01;
NUMBER_OF_DOMAINS(iproj) = ndom	;
%---------------------------------------
iproj = 3;
NAMEPROJECTS{iproj} = [NAME_PROJECTS,'_',NAMES_PROJECTS_PP{iproj}] ;
IMPOSED_DISP.RIGHT_END(iproj,1)= 0.01;
NUMBER_OF_DOMAINS(iproj) = ndom	;
%---------------------------------------
iproj = 4;
NAMEPROJECTS{iproj} = [NAME_PROJECTS,'_',NAMES_PROJECTS_PP{iproj}] ;
IMPOSED_DISP.RIGHT_END(iproj,4)= 0.01;
NUMBER_OF_DOMAINS(iproj) = ndom	;
%---------------------------------------
iproj = 5;
NAMEPROJECTS{iproj} = [NAME_PROJECTS,'_',NAMES_PROJECTS_PP{iproj}] ;
IMPOSED_DISP.RIGHT_END(iproj,5)= 0.01;
NUMBER_OF_DOMAINS(iproj) = ndom	;
%---------------------------------------
iproj = 6;
NAMEPROJECTS{iproj} = [NAME_PROJECTS,'_',NAMES_PROJECTS_PP{iproj}] ;
IMPOSED_DISP.RIGHT_END(iproj,6)= 0.01;
NUMBER_OF_DOMAINS(iproj) = ndom	;
%---------------------------------------


%

% DILATATION_SLICES_x = linspace(1,1,NUMBER_OF_DOMAINS) ; % Dilatation of slices in the x-direction
%DILATATION_SLICES_x = ones(1,NUMBER_OF_DOMAINS);


DATA_TYPESOLVER = 0; % Conjugated gradient if = 1.
DATA_niterCONJG = 100000 ; % Number of iterations  conjugated gradient


FOLDER_SAVE = [FOLDER,'/DATAFE_TRAIN/'] ;
if exist(FOLDER_SAVE)==0;     mkdir(FOLDER_SAVE) ; end

DATA.RECALCULATE_STIFFNESS = 1 ;

if ~exist('FE_ELASTOSTATIC','file')
    addpath([EXECUTABLE_FOLDER,'FE_CODE']) ;
end
PROJECT_REFERENCE = [] ;
FLUCTUATIONS = [] ;
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
        
        eval(NAME_DATA_REF) ;
        
        cd(EXECUTABLE_FOLDER) ;
        DATA.INPUTDATAfile = [FOLDER,'/',NAMEPROJECTS{iprojects}] ;
        FUNinput.INPUTS.DISPLACEMENT_COND_FILE =TYPE_BOUNDARY_CONDITIONS ;
        FUNinput.INPUTS.BEAMLOADS.LOAD_UNIFORM = cell(1,5) ;
        FUNinput.INPUTS.BEAMLOADS.LOAD_UNIFORM(:) = {zeros(NUMBER_OF_DOMAINS(iprojects),3) };
        % FE CODE
        DATA.MakeMeshByRepetition.nDOM = NUMBER_OF_DOMAINS(iprojects) ;
        %   DATA.MakeMeshByRepetition.DILATATION_SLICES_x = DILATATION_SLICES_x ;
        MATERIAL.PLY(1).E = E ;
        MATERIAL.PLY(1).nu = nu ;
        MATERIAL.PLY(1).G = G ;
        MATERIAL.PLY(1).typePROBLEM = typePROBLEM ;
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
        
        
        if iprojects == 5  
            
            DATA.FLUCTUATION_MODES_BENDING = FLUCTUATIONS(:,1)  ;
        elseif iprojects == 6 
            DATA.FLUCTUATION_MODES_BENDING = FLUCTUATIONS(:,2)  ;
        else
            DATA.FLUCTUATION_MODES_BENDING  = [] ;
        end
        
        DATAOUT = FE_ELASTOSTATIC(FUNinput,DATA) ;
        
        
        DATA_FE_WS{iprojects} =  DATAOUT.nameWORKSPACE ;
        
        if  iprojects == 1 || iprojects == 2
            disp('Computing fluctuations')
            disp('_______________________________')
            
            NODES_A = DATAOUT.DOMAINVAR.NODES_faces12{1} ;
            ndim = 3;
            DOFsA = small2large(NODES_A,ndim) ;
            R = DATAOUT.DOMAINVAR.RigidBodyModes ;
            a_A = (FUNinput.INPUTS.DISPLACEMENTS_ENDS_RIGID_BODY{1}) ;
            
            U_A = DATAOUT.d(DOFsA) - R*a_A' ;
            TOL = 1e-6 ;
            if norm(U_A)/norm(DATAOUT.d(DOFsA))<TOL
                U_A =[]  ; % Negligible fluctuations
            end
            FLUCTUATIONS  = [FLUCTUATIONS U_A ];
            
            
            if iprojects == 2 & ~isempty(FLUCTUATIONS)
                load(DATAOUT.nameWORKSPACE,'COOR','TypeElementB','CONNECTb','posgp')
                CNref =  RenumberConnectivities(CONNECTb{1},1:length(NODES_A)) ;
                COOR = COOR(NODES_A,:) ;
                
                DATAMODES = [] ; DOFl = [] ;
                NAME_MODES = 'fluct_pure_b_yz' ;
                LEGENDG= ['fluct_pure_b_yz'] ;
                posgp = [] ;
                NAME_MODES = [DATAOUT.nameWORKSPACE(1:end-4),LEGENDG ];
                
                GidPostProcessModes_dom(COOR,CNref,TypeElementB,FLUCTUATIONS,posgp,NAME_MODES,DATAMODES,LEGENDG);
                % DATA.FLUCTUATION_MODES_BENDING = FLUCTUATIONS ;
                %  FLUCTUATIONS =
                DATA.FLUCTUATION_MODES_BENDING = FLUCTUATIONS ; 
            end
            
        else
            % Checking that the implementation is correct 
            load(DATAOUT.nameWORKSPACE,'COOR','TypeElementB','CONNECTb','posgp')
            [CentroidFA,AREA,Mst] =...
                CentroidGeometricMassMatrixNEW(COOR,NODES_A,CONNECTb{1},TypeElementB) ;
            M = sparse(3*size(Mst)) ;
            for idim =1:3
                INDLOC =idim:3:size(R,1) ;
                M(INDLOC,INDLOC) =  Mst ;
            end
             NODES_B = DATAOUT.DOMAINVAR.NODES_faces12{2} ;
             DOFsB = small2large(NODES_B,ndim) ;
             dA = DATAOUT.d(DOFsA) ; 
             dB = DATAOUT.d(DOFsB) ; 
               a_A = (FUNinput.INPUTS.DISPLACEMENTS_ENDS_RIGID_BODY{1}) ;
                 a_B = (FUNinput.INPUTS.DISPLACEMENTS_ENDS_RIGID_BODY{2}) ;
                 Ub  = FLUCTUATIONS ; 
                 dCb  = (Ub'*M*Ub)\(Ub'*M*(dA-dB)) ; 
                 
                 Q = Ub*((Ub'*M*Ub)\(Ub'*M)) ; 
                 
                 compr = (dA-dB) - Q*(dA-dB) ; 
            
            
        end
        
        
        
        
        %
        cd(FOLDER) ;
        %    DATA.RECALCULATE_STIFFNESS = 0 ;
        
        
    end
end

save([FOLDER_SAVE,filesep,[NAMEFILE,'_FEWS.mat']],'DATA_FE_WS')
