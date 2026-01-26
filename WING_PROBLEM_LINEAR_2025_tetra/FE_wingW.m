clc
clear all

if ~exist('EIFEvectorizedGEN1')
    addpath(genpath(getenv('MATLAB_CODES_FEHROM')))  ; %
end

% LINEAR REGIME
NAMEMESH =  'FEwingW' ;
DATALOC.NameFileMeshDATA =[cd,filesep,'GIDPRE/',NAMEMESH,'.msh'] ;
BOUNDARY_CONDIDITONS = 'BCfe_CWeight' ;
MATERIAL_FILE  = 'MaterialDATA_elas3D_2mat' ;
%DATALOC.OnlyComputeMassAndStiffnessMatricesSmallStrains = 1;
DATALOC1.FACE_TO_ANALIZE = 1 ;
DATALOC1.DOF_TO_ANALIZE = 2 ;

COMPUTE_FE   =0;

DATALOC1.LEGEND_LOC = 'FEM' ;

if COMPUTE_FE == 1
    FEvectorized_GEN1('INPUTS_FE_LOC',BOUNDARY_CONDIDITONS,MATERIAL_FILE,DATALOC) ;
    DATALOC = [] ;
    ExtractINFOplotNL_surfUD(BOUNDARY_CONDIDITONS,DATALOC1) ; %
else
    % Reactions
    
    ExtractINFOplotNL_surfUD(BOUNDARY_CONDIDITONS,DATALOC1) ;
    
    % VIBRATION MODES
end
