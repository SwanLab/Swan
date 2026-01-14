clc
clear all

if ~exist('Quadratic3N')
    addpath(genpath('CODE'))  ; %
end
% LINEAR REGIME
NAMEMESH = 'EIFE_wingC100';   %'Beam02x2_5_10x100q' ;
DATALOC.NameFileMeshDATA =[cd,filesep,'GIDPRE/',NAMEMESH,'.msh'] ;
%load('ListElements.m') ; ListElements = ListElements(:,1) ;
DATALOC.ElementToStudy =[1:6,8,10,40:41,90:100] ;  ; % [14,26];
BOUNDARY_CONDIDITONS = 'BCeifem_CWeight' ;
MATERIAL_FILE  = 'MATERIAL_EIFE_2matEL' ;
DATALOC.FineScaleDomainsToShow =     DATALOC.ElementToStudy  ; %
NAMEWS_save_strainHISTORY = 'CSTRAIN_histories/Element14.mat' ;
DATALOC.ONLY_PRINT_GID  =0; % Only post-process

%DATALOC.OnlyComputeMassAndStiffnessMatricesSmallStrains = 1;
DATALOC1.FACE_TO_ANALIZE = 1 ;
DATALOC1.DOF_TO_ANALIZE = 2 ;
DATALOC1.LEGEND_LOC = 'EIFEM' ;

COMPUTE_EIFEM   =1;


if COMPUTE_EIFEM == 1
    EIFEvectorizedGEN1('INPUTS_EIFE_shapeBUBloc',BOUNDARY_CONDIDITONS,MATERIAL_FILE,DATALOC) ;
    
    ExtractINFOplotNL_surfUD(BOUNDARY_CONDIDITONS,DATALOC1) ;
else
    % Reactions
    ExtractINFOplotNL_surfUD(BOUNDARY_CONDIDITONS,DATALOC1) ;
end
