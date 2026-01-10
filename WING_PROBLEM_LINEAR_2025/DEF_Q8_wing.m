function DATAOUT=    DEF_Q8_wing(DATAoffline)
if nargin == 0
    DATAoffline = [] ; 
end
% ---------------------------------------------------------------------
% Input parameter specific for a given project
% Training Q8 plate elements
% ---------------------------------------------------------------------------------------------

% ---------------------------------------------------------------------
% 0) TEMPLATE INPUT function
% This is the actual input data file. The outputs of the present input
% file(function ) are fed into this function in order to produce data that
% can be processed by the FE software
DATAOUT.InputDataFile_FE = 'INPUTS_TRAIN_Q8plate' ;

%1) Label identifying parameter study
% --------------------------------------------
DATAOUT.NameParamStudyLOC = ['wingQ8'] ;

%2) Order interpolation macro-element  (this is the "TRAINING" domain, for instance, a 3x3 cell)
% --------------------------------------------
NAME_MACRO_ELEMENT_MESH = [] ;%  Not necessary for this type of training 
DATAOUT.COARSE_ELEMENT_FESHAPE.NAMEFILEMESH = [] ; 
DATAOUT.TypeProblem = '3D' ; 

% DATAOUT.COARSE_ELEMENT_FESHAPE.NAMEFILEMESH =[cd,filesep,'GIDPRE/',NAME_MACRO_ELEMENT_MESH,'.msh'] ;
% DATALOC_trainingDOM.NameFileMeshDATA =  DATAOUT.COARSE_ELEMENT_FESHAPE.NAMEFILEMESH  ;
% DATALOC_trainingDOM.READ_MATERIAL_COLUMN = 0 ;
% DATALOC_trainingDOM.RenumberElementsForEficiency = 0 ;
% DATAOUT.MESHcoarse = GeometryMesh(DATALOC_trainingDOM) ;


%3 ) TIME STEPS  (Dummy in elastic problems )
% ------------------------%
DATAOUT.t0 = 0 ; DATAOUT.tEND = 1 ; ntimes = 2 ;
DATAOUT.DATA_STEPS = linspace(DATAOUT.t0,DATAOUT.tEND,ntimes) ;
% For printing in GID
DATAOUT.STEPS_print_FREQ = 1 ; % Only
% Steps to print in GID (SUBSET OF DATA.STEPS)
DATAOUT.DATA_PRINT_NSTEPS =  unique([1:DATAOUT.STEPS_print_FREQ:length(DATAOUT.DATA_STEPS),length(DATAOUT.DATA_STEPS)]) ;


% 4) FINITE ELEMENT MESH
% ------------------------------
NAMELOCMESH ='NACA2s_trainingdom';
DATAOUT.NameFileMeshDATA =[cd,filesep,'GIDPRE/',NAMELOCMESH,'.msh'] ;
% FE mesh of the studied geometry. Remember to load
% problem type: /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/PROBLEMTYPES_GID/PROBLEM_TYPE_SIMPLE.gid
% To generate the files needed by matlab,
% remember to export both the mesh and the conditions data. (Files >
% Export > Gid Mesh) and (Files > Export > Calculation File)
% Furthermore, you must identify the domain under study, the 6 surfaces
% defining the training domain, and the 6 surfaces defining the trained
% domain itself
% SURF1 - XMIN, SURF2 -XMAX, .... SURF6-ZMAX (TRAINING DOMAIN)
% SURF7- XMIN, ....SURF12 .... ZMAX  (TRAINED DOMAIN)

% Mesh parent domain 
MESH_PARENT_DOMAIN_COARSE ='NACA2s_Q8_trainedDOMAIN';
DATAOUT.MESH_PARENT_DOMAIN_COARSE =[cd,filesep,'GIDPRE/',MESH_PARENT_DOMAIN_COARSE,'.msh'] ;

 


% 3) PARAMETERIZATION DIRICHLET BOUNDARY CONDITION/NEUMANN CONDITIONS (PARAMETER SPACE)
% --------------------------------------------------------------------
% TRAINING WITH 48 TESTS 
 
AMPLITUDE_alltests = 0.1 ;
index_SURFACES_PRESCRIBED_DISPLACEMENTS =[1,2] ;
%[DATAOUT] = Q8_for_beam_12testsF(AMPLITUDE_alltests,DATAOUT,index_SURFACES_PRESCRIBED_DISPLACEMENTS) ; 

%[DATAOUT] = Q8_for_beam_18tests(AMPLITUDE_alltests,DATAOUT,index_SURFACES_PRESCRIBED_DISPLACEMENTS)
 
 %[DATAOUT] = Q8_for_beam_8tests(AMPLITUDE_alltests,DATAOUT,index_SURFACES_PRESCRIBED_DISPLACEMENTS)  ; 
 
  [DATAOUT] = Q8_for_beam_6tests(AMPLITUDE_alltests,DATAOUT,index_SURFACES_PRESCRIBED_DISPLACEMENTS)  ; 

 
%
DATAOUT.LabelEntitiesDefiningBoundary = [1,2,3,4,5,6] ; % Label lines defining Dirichlet boundaries
% -------------------

DATAOUT.TypeFunctionDisplacementInterfaces =  'HEXAHEDRA_LINEAR' ;
DATAOUT.INCLUDE_FLUCTUATION_DOFS_IN_COARSE_SCALE = 0;  % = 0 , IT IGNORES FLUCTUATIONS


% 4.1) MATERIAL
DATAOUT.FILE_MATERIAL_DATA = 'MaterialDATA_elas3D_4mat' ;



% 6)   KINEMATICS
% -------------------
DATAOUT.SMALL_STRAIN_KINEMATICS= 1;



ENRICHED_GAUSS_RULE =0;
if ENRICHED_GAUSS_RULE == 1
    ngausLOC = [4,4];
    [~, ~, xGAUSS, weights]  =TensorProd2Ddiscr(ngausLOC) ;
    [xx,yy]  = meshgrid(xGAUSS{1},xGAUSS{2}) ;
    xx = xx(:) ; yy = yy(:);
    posgp= [xx'; yy'] ;
    DATAOUT.posgp_given  = posgp ;
    DATAOUT.weights_given  = weights ;
end






