function  ExtractModesPartition2018(DATA_TRAINING,DATAIN)
% Modal analysis domain-wise,
%dbstop('4')
if nargin ==0
    load('tmp.mat')
end
if exist('ExtractDisplMatrix')~=2
    addpath('FE_CODE') ;
end
if exist('BeamStiffnessMatrix')~=2
    addpath('FE_CODE/BeamROM') ;
end
if exist('SVDT')==0
    addpath('SVDlibrary')
end
if exist('EmpiricalCubatureMethod')==0
    addpath(genpath('SVDlibrary'))
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT inputs
DATAIN= DefaultField(DATAIN,'NMODES_SHOW',[]) ;
DATAIN = DefaultField(DATAIN,'CUBATURE',[]) ;

DATAIN = DefaultField(DATAIN,'CUBATURE',[]) ;
DATAIN.CUBATURE = DefaultField(DATAIN.CUBATURE,'CUBATURE',0) ; % Efficient integration scheme
DATAIN = DefaultField(DATAIN,'TOLERANCE_SVD_DISPLACEMENTS',[]) ;  % If not-empty, the number of modes are determined by this tolerance
%%%%%%%%%%%%%%%%%%%%%%%%% END INPUTS %%%%%%%%%%%%%%%%%%%
MSG = {} ; 
% -------------------------------------
% CONSTRUCTING BASIS MATRICES
% ------------------------------------
DATAIN = DefaultField(DATAIN,'COMPUTE_MODES_AGAIN',1) ;
if DATAIN.COMPUTE_MODES_AGAIN == 1
    % ---------------------------------------------------------------------
    DATAIN = DefaultField(DATAIN,'ISNONLINEAR',0)  ;
    %    if   DATAIN.ISNONLINEAR == 0
    % From 28-May-2019, the same function is used for both cases
    [BASES,DATA_REFMESH,nBASES_BEAM,DATAsnap,MSG ]= BasisU_R_def_computation(DATA_TRAINING,DATAIN,MSG) ;
    %   else
    %     [BASES,DATA_REFMESH,nBASES_BEAM,DATAsnap ]= BasisU_R_def_NONLINEAR(DATA_TRAINING,DATAIN) ;
    %  end
    % -------------------------------------------------------------------------
    
    save(DATAIN.NAME_WS_MODES,'BASES','DATA_REFMESH','nBASES_BEAM','DATAsnap','-append');
else
    load(DATAIN.NAME_WS_MODES,'BASES','DATA_REFMESH','nBASES_BEAM','DATAsnap')
end

DATAIN = DefaultField(DATAIN,'APPROACH_STRAINING_MODES_GREATER_THAN_REACTION_MODES',0)  ; 

if DATAIN.APPROACH_STRAINING_MODES_GREATER_THAN_REACTION_MODES == 0
    % Before May-13th-2020  (including Paper's result)
    [DATAROM,MSG] = BeamStiffnessMatrix(BASES,DATA_REFMESH,DATAIN,nBASES_BEAM,DATAsnap,MSG) ;
else
%    error('This option should be explored more in depth... (13thMAy-2020)')
    % After May-13th-2020. GOAL: Explore what happens if the number of
    % straining modes are greater than the number of reaction modes
    [DATAROM,MSG] = BeamStiffnessMatrix_MODESnu_g_nr(BASES,DATA_REFMESH,DATAIN,nBASES_BEAM,DATAsnap,MSG) ;
end




save(DATAIN.NAME_WS_MODES,'DATAROM','-append');


%%% COMPUTING REDUCED (BEAM STIFFNESS MATRIX) OF THE REFERENCE ELEMENT

DATAIN = DefaultField(DATAIN,'COMPUTE_STRESSES_AT_REDUCED_POINTS',1) ;

%if  DATAIN.COMPUTE_STRESSES_AT_REDUCED_POINTS == 1
% Empirical Cubature Method ---determination of integration points

DATAIN.CUBATURE = DefaultField(DATAIN.CUBATURE,'CONTINUOUS_ECM',0) ;  % Continuous ECM 
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/ContinuousEmpiricalCubatureM/Paper_hernandez2021ecm/09_ASSESSMENTpaper/MultiscaleHROM/READMEmulti.mlx

if DATAIN.CUBATURE.CONTINUOUS_ECM == 0
[HROMVAR,MSG ]= ReducedSetIntegrationPoints(DATAIN,DATAROM,BASES,DATA_REFMESH,MSG) ;
else
% 3th November 2022
[HROMVAR,MSG ]= ReducedSetIntegrationPointsCECM(DATAIN,DATAROM,BASES,DATA_REFMESH,MSG) ;
error('The remaining part of the code has to be adapted for this option (Nov-2022)')
end

% Reduced-order operators required for nonlinear analysis
HROMVAR = ROMoperatNONL(DATAROM,BASES,DATA_REFMESH,DATAIN,HROMVAR) ;


DATAROM.HROMVAR = HROMVAR ;
save(DATAIN.NAME_WS_MODES,'DATAROM','DATAIN','-append')
%end

% Nonlinear

% Saving information 
% -----------------------
[LOCAL_FOLDER,NAMEFILELOC, EXTENSION]= fileparts(DATAIN.NAME_WS_MODES) ; 
 DATAIN = DefaultField(DATAIN,'LOCAL_FOLDER',LOCAL_FOLDER) ;
 LOCAL_FOLDER = DATAIN.LOCAL_FOLDER; 
 REPORTS_FOLDER = [LOCAL_FOLDER,filesep,'REPORTS',filesep] ; 
if  ~exist(REPORTS_FOLDER,'file')
    mkdir(REPORTS_FOLDER) ; 
end

aaa = datestr(datetime) ;

FILE_BATCH  = [REPORTS_FOLDER,DATAIN.NAME_project,'_',aaa,'.txt'] ;
fid =fopen(FILE_BATCH,'w');
for i = 1:length(MSG)
    for j = 1:size(MSG{i},1)
        fprintf(fid,[MSG{i}(j,:),'\n']);
    end
end
fod =fclose(fid);

open(FILE_BATCH);



FILE_COPY = [REPORTS_FOLDER,DATAIN.NAME_project,'_',aaa,'.m'] ;
copyfile([LOCAL_FOLDER,filesep,DATAIN.NAME_project,'.m'],FILE_COPY)  ; 

