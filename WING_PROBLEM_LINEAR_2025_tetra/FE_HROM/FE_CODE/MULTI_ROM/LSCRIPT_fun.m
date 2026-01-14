function LSCRIPT_fun(EXECUTABLE_FOLDER,NAME_GEOMETRY_MATERIAL_DATA,NAME_LOAD_DATA,RUN_MODES,RUN_FE)

if nargin == 0
    % Input data
    EXECUTABLE_FOLDER = '/home/joaquin/Desktop/CURRENT_TASKS/COMPOSITE_MATERIALS_DOCENCIA/APUNTES_LATEX/DOCUMENTOS_anexos/MATLAB/ELASTOSTATIC_GEN/' ;
    NAME_GEOMETRY_MATERIAL_DATA ='DATA_INPUT_e11' ;  'DATA_INPUT_e10' ;
    NAME_LOAD_DATA =  'LOADS_biempotrada' ; 'LOADS_torsion' ;   ;'LOADS_empotrada' ;
    
    
    
    
end

% ONLINE OPERATIONS
% --------------------------------------
RUN_ROM =0;  % Run online operations

% SLICES
% ----------------------------
RUN_FE.SLICES.BEAM =1;  % Run FE analyses for slices (``beam" training)
% Determination of modes -- BEAMS
% --------------------------------
RUN_MODES.SLICES.BEAM =1;        % Compute reduced-order operators
COMPUTE_MODES_AGAIN.SLICES.BEAM =1; % Compute modes. If = 0, it means they have been already calculated

% END INPUT DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NO_OFFLINE_OPERATIONS = 0;
% -------------------------------
% JOINTS
% --------------------------------
RUN_FE.JOINTS = [0 0];  % Run FE analyses for joints
% Determination of modes
RUN_MODES.JOINTS = [0 0];
COMPUTE_MODES_AGAIN.JOINTS= [0 0];


FOLDER = cd;
MFILENAME = NAME_GEOMETRY_MATERIAL_DATA;

if NO_OFFLINE_OPERATIONS == 1
    GENERATE_JOINTS = 0 ;
    
    for ibeams = 1:length( RUN_FE.SLICES.BEAM )
        RUN_FE.SLICES.BEAM(ibeams) = 0 ;
    end
    for ibeams = 1:length( RUN_FE.JOINTS)
        RUN_FE.JOINTS(ibeams) = 0 ;
    end
    
    for ibeams = 1:length( RUN_MODES.SLICES.BEAM )
        RUN_MODES.SLICES.BEAM(ibeams) = 0 ;
    end
    
    for ibeams = 1:length( COMPUTE_MODES_AGAIN.SLICES.BEAM )
        COMPUTE_MODES_AGAIN.SLICES.BEAM(ibeams) = 0 ;
    end
    
    for ibeams = 1:length( RUN_MODES.JOINTS )
        RUN_MODES.JOINTS(ibeams) = 0 ;
    end
    
    for ibeams = 1:length( COMPUTE_MODES_AGAIN.JOINTS )
        COMPUTE_MODES_AGAIN.JOINTS(ibeams) = 0 ;
    end
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DATARUN.OPEN_GID_ROMRES = 1;
PLAY_SOUND = 0;


GENERATE_JOINTS =1;   %  Guide the user in creating  joints
PLOT_3D_STRUCTURE_CHECK =0;  % PLot the entire 3D structure (for checking purposes)
if ~exist('INPUTS_GENERAL','file') ; addpath([EXECUTABLE_FOLDER,'/FE_CODE/']);  end
if ~exist('SVD_dom','file') ; addpath([EXECUTABLE_FOLDER,'/FE_CODE/MULTILEARN/']);  end
if ~exist('GeometryStructure','file') ; addpath([EXECUTABLE_FOLDER,'/FE_CODE/BeamROM/']);  end
Training_and_ROM_straightbeam(NAME_GEOMETRY_MATERIAL_DATA,FOLDER,EXECUTABLE_FOLDER,MFILENAME,...
    RUN_FE,GENERATE_JOINTS,RUN_MODES,RUN_ROM,COMPUTE_MODES_AGAIN,...
    NAME_LOAD_DATA,PLOT_3D_STRUCTURE_CHECK,DATARUN) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if PLAY_SOUND ==1
    load handel.mat;
    sound(y, Fs);
end
cd(FOLDER)