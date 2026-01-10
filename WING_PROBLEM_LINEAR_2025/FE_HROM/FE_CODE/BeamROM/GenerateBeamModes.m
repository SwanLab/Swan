function GenerateBeamModes(NAME_MODES,COMPUTE_MODES_AGAIN,TOL_SINGULAR_VALUES_Hqr,...
    NAME_ROOT_FE,EXECUTABLE_FOLDER,DATAIN,NonBeamProjects,DATARUN,BeamProjects)

if nargin == 0
    load('tmp2.mat')
elseif nargin == 5
    DATAIN = [] ;
end

DATAIN = DefaultField(DATAIN,'SLICES_SEL',[]) ;
DATAIN.TOL_SINGULAR_VALUES_Hqr = TOL_SINGULAR_VALUES_Hqr ;
FOLDER = [cd] ;
DATAIN.PLOT_RIGHT_SINGULAR_VECTORS = 10 ;
 
% Names FE input files    % Number of disp. ,modes   % Slices to be sampled
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NAME_PROJ_LOC = {'axial','torsion','pbend_y','pbend_z','sbend_y','sbend_z'} ;
if DATARUN.ndim == 3
    NAME_PROJ_LOC = {'pbend_y','pbend_z','axial','torsion','sbend_y','sbend_z'} ;
    
    [~,NAME_PROJ_LOC,~,~,~] ...
        = NameProjectsBeam([],[],[],DATARUN) ;
    
    
else
    %   NAME_PROJ_LOC = {'axial','pbend_z','sbend_z'} ;
    [~,NAME_PROJ_LOC,~,~,~] ...
        = NameProjectsBeam2D([],[],[],DATARUN) ;
end
DATAIN = DefaultField(DATAIN,'NMODES_PROJECT_DISPLACEMENTS',cell(size(NAME_PROJ_LOC)))  ;
DATAIN = DefaultField(DATAIN,'NMODES_PROJECT_REACTIONS',cell(size(NAME_PROJ_LOC)))  ;
DATAIN = DefaultField(DATAIN,'NMODES_PROJECT_STRESSES',cell(size(NAME_PROJ_LOC)))  ;


%DATAIN = DefaultField(DATAIN,'UnifiedApproachForModes',0) ; %  New approach, 10-dec-2018 (No distinction between beam modes and non-beam modes)

if  DATAIN.UnifiedApproachForModes == 1
    %  New approach, 10-dec-2018 (No distinction between beam modes and non-beam modes)
    DATAIN = DefaultField(DATAIN,'TOL_PROJECT_SVD',[]) ;
    DATAIN.TOL_PROJECT_SVD = DefaultField(DATAIN.TOL_PROJECT_SVD,'DISP',1e-6) ;
    DATAIN.TOL_PROJECT_SVD = DefaultField(DATAIN.TOL_PROJECT_SVD,'REACT',1e-6) ;
    DATAIN.TOL_PROJECT_SVD = DefaultField(DATAIN.TOL_PROJECT_SVD,'STRESS',1e-6) ;
    TOL_PROJECT_SVD = DATAIN.TOL_PROJECT_SVD ;
    
else
    % No longer used 
    error('Option no longer used (Aug-2019)')
    TOL_PROJECT_SVD.DISP = [] ;
    TOL_PROJECT_SVD.REACT = [] ;
    TOL_PROJECT_SVD.STRESS = [] ;
end

SLICES_SEL = DATAIN.SLICES_SEL ;
IS_BEAM_MODE = ones(size(NAME_PROJ_LOC)) ;
for iproj = 1:length(NAME_PROJ_LOC)
    NAMEPROJECTS{iproj} = [NAME_ROOT_FE,NAME_PROJ_LOC{iproj}];
    nU{iproj}= DATAIN.NMODES_PROJECT_DISPLACEMENTS{iproj};
    nR{iproj}= DATAIN.NMODES_PROJECT_REACTIONS{iproj};
    nS{iproj}= DATAIN.NMODES_PROJECT_STRESSES{iproj};;
    nINTF{iproj}= [];
    TOLERANCE_SVD_SLICE{iproj}= TOL_PROJECT_SVD.DISP;
    TOLERANCE_SVD_SLICE_REACT{iproj}= TOL_PROJECT_SVD.REACT;
    TOLERANCE_SVD_SLICE_STRESS{iproj}= TOL_PROJECT_SVD.STRESS;
    SLICES{iproj} =SLICES_SEL ;
    TIME_STEPS{iproj} = []  ;
end

% iproj = 7;
% NAMEPROJECTS{iproj} = 'TEST_pressure' ;
% nU{iproj}= [10];
% nR{iproj}= [];
% SLICES{iproj} =[] ;

DATAIN.BeamProjects = BeamProjects ;


%IS_BEAM_MODE(7) = 0 ;


%%%%%%%%%%%%%%%%%


TOTAL_NUMBER_OF_MODES = [] ;

%

%nR = nU ;  % Number of reaction modes
% If SLICE{i} is empty, then all slices are included in the SVD
if ~isempty(NAME_PROJ_LOC)
    DATA_TRAINING = cell(size(NAMEPROJECTS)) ;  % % Full Path Names .mat files containing FE data (Stiffness matrix...)
    for itraj = 1:length(DATA_TRAINING)
        %  These paths are specified when running the FE analyses (FE_TRAIN_fast.m, in this case)
        DATA_TRAINING{itraj} = [ FOLDER,filesep,'DATAFE_TRAIN',filesep,NAMEPROJECTS{itraj},'.mat'] ; %
    end
    
else
    DATA_TRAINING = [] ;
    nU = [] ; %DATAIN.NMODES_PROJECT_DISPLACEMENTS{iproj};
    nR  = [] ;%  DATAIN.NMODES_PROJECT_REACTIONS{iproj};
    nS =[] ; %DATAIN.NMODES_PROJECT_STRESSES{iproj};;
    nINTF= [];
    TOLERANCE_SVD_SLICE =[] ; %[ ; %{iproj}= TOL_PROJECT_SVD.DISP;
    TOLERANCE_SVD_SLICE_REACT=  [] ; %TOL_PROJECT_SVD.REACT;
    TOLERANCE_SVD_SLICE_STRESS= [] ; % TOL_PROJECT_SVD.STRESS;
    SLICES =[] ;
    TIME_STEPS = [] ;
end




if ~isempty(NonBeamProjects)
    
    NAMEWSloc =    NonBeamProjects.NAMEWS  ;   % Absolute path
    nPROJADD = length(NAMEWSloc)
    NonBeamProjects  =DefaultField(NonBeamProjects,'nU',cell(nPROJADD,1)) ;
    NonBeamProjects  =DefaultField(NonBeamProjects,'nR',cell(nPROJADD,1)) ;
    NonBeamProjects  =DefaultField(NonBeamProjects,'nS',cell(nPROJADD,1)) ;
    NonBeamProjects  =DefaultField(NonBeamProjects,'nINTF',cell(nPROJADD,1)) ;
    NonBeamProjects  =DefaultField(NonBeamProjects,'TOLERANCE_SVD_SLICE',cell(nPROJADD,1)) ;
    NonBeamProjects  =DefaultField(NonBeamProjects,'TOLERANCE_SVD_SLICE_STRESS',cell(nPROJADD,1)) ;
    NonBeamProjects  =DefaultField(NonBeamProjects,'TOLERANCE_SVD_SLICE_REACT',cell(nPROJADD,1)) ;
    
    
    NonBeamProjects  =DefaultField(NonBeamProjects,'SLICES',cell(nPROJADD,1)) ;
    NonBeamProjects  =DefaultField(NonBeamProjects,'TIME_STEPS',cell(nPROJADD,1)) ;
    
    for iaddproj = 1:length(NAMEWSloc)
        DATA_TRAINING{end+1} = NAMEWSloc{iaddproj} ;
        nU{end+1} = NonBeamProjects.nU{iaddproj}  ;
        nR{end+1} = NonBeamProjects.nR{iaddproj}  ;
        nS{end+1} = NonBeamProjects.nS{iaddproj}  ;
        nINTF{end+1} = NonBeamProjects.nINTF{iaddproj}  ;
        SLICES{end+1} =NonBeamProjects.SLICES{iaddproj}  ;
        TIME_STEPS{end+1} =NonBeamProjects.TIME_STEPS{iaddproj}  ;
        TOLERANCE_SVD_SLICE{end+1} =NonBeamProjects.TOLERANCE_SVD_SLICE{iaddproj}  ;
        TOLERANCE_SVD_SLICE_REACT{end+1} =NonBeamProjects.TOLERANCE_SVD_SLICE_REACT{iaddproj}  ;
        TOLERANCE_SVD_SLICE_STRESS{end+1} =NonBeamProjects.TOLERANCE_SVD_SLICE_STRESS{iaddproj}  ;
        if DATAIN.UnifiedApproachForModes == 0
            IS_BEAM_MODE(end+1) = 0 ;
        else
            IS_BEAM_MODE(end+1) = 1 ;
        end
    end
    
    
    
    
end




% % NAME BINARY FILE (.mat) WHERE MODES  TO BE COMPUTED WILL BE STORED
if ~exist('MODES','dir')
    mkdir('MODES') ;
end
if isempty(NAME_MODES)
    NAMEFILE = mfilename ;
    
else
    NAMEFILE = NAME_MODES ;
end
DATAIN.NAME_WS_MODES = [cd,filesep,'MODES',filesep,'MODES_',NAMEFILE,'.mat'];

%DATAIN.NAMEWS = [cd,'/','MODES_',NAMEFILE,'.mat'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of modes of each project for deformational displacements,
DATAIN.NMODES_PROJECT_DISPLACEMENTS =  nU;
DATAIN.NMODES_PROJECT_STRESSES =  nS;
DATAIN.TOLERANCE_SVD_SLICE = TOLERANCE_SVD_SLICE ;
DATAIN.TOLERANCE_SVD_SLICE_REACT = TOLERANCE_SVD_SLICE_REACT ;
DATAIN.TOLERANCE_SVD_SLICE_STRESS = TOLERANCE_SVD_SLICE_STRESS ;
% Number of reaction modes
DATAIN.NMODES_PROJECT_REACTIONS =   nR;
DATAIN.NMODES_PROJECT_INTERFACE =   nINTF;

DATAIN.DOMAINS_TO_INCLUDE_TRAINING = SLICES;
DATAIN.TIME_STEPS_TO_INCLUDE_TRAINING = TIME_STEPS;
DATAIN.NMODES_TRUNCATE = TOTAL_NUMBER_OF_MODES ;
DATAIN.IS_BEAM_MODE = IS_BEAM_MODE ;
DATAIN.COMPUTE_MODES_AGAIN = COMPUTE_MODES_AGAIN ;

%%% Information displacement face B in axial, shear.... tests (for beam modes)
DATAIN.NameWs_displacementB_beams =    [ FOLDER,filesep,'DATAFE_TRAIN',filesep,NAME_ROOT_FE,'global.mat'] ; % [NAME_ROOT_FE,'global.mat'] ;

%%% END INPUT DATA %----------------------------------
cd(EXECUTABLE_FOLDER)
ExtractModesPartition2018(DATA_TRAINING,DATAIN) ;
cd(FOLDER)