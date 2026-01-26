clc
clear all
% Finite element code for composite materials
%  Technical University of Catalonia, DIC-2015
%  % Joaquín A. HERNANDEZ, jhortega@cimne.upc.edu
% --------------------------------------------------------------------------------------------

NameFileMesh ='wing1_1000el.msh'; 'lam10lay_8000.msh' ; 'shear.msh' ; 'wing1.msh' ;'wing3.msh' ;'lam10lay_8000.msh' ;  'wing1_1000el.msh'; 'wing1_1000el.msh';'wing1.msh' ;  'lam10lay_8000.msh' ;   % Geometry of the laminate
%NameWS_unitcell = 'DATAWS/Celas_iSO.msh.mat' ;'DATAWS/Celas_ejem1.msh.mat' ;   % 
NameWS_unitcell ='DATAWS/Celas_ejem1.msh.mat' ;
NameWS_failureE = 'DATAWS/FAILURE_ejem1.msh.mat' ;
DATA.NUMBER_ELEMENTS_FAIL = 2;
DATA.PLOT_FAILURE_POINTS = 1 ; 
DATA.CALCULATE_MASSMATRIX = 1 ; 
DATA.STORE_STIFFNESS = 1 ; %Store stiffness matrix

% Elasticity matrix for each  layer (and angle of rotation)
% The name of the .mat containg such information is to be provided  .
% -----------------------
%
nplies = 8 ; 
ANG_FIB =    [0 90 0 0 0 0 90 0]; zeros(nplies,1) ; 
for iply = 1:length(ANG_FIB)
    MATERIAL.PLY(iply).NAMEWS = NameWS_unitcell;  % Store the Celas in this folder
    MATERIAL.PLY(iply).ANGLE = ANG_FIB(iply) ;%;  % Angle subtended by the fiber and the x-axis
    MATERIAL.PLY(iply).NAMEWS_FAIL = NameWS_failureE ;
end
 
DATA.TYPESOLVER = 1; % Conjugated gradient
DATA.niterCONJG = 3000 ; % Number of iterations  conjugated gradient
% INPUT DATA
% ----------
FUNinput.NAME = 'INPUT_LAMINATES_GEN' ; % NAME OF THE FUNCTION THAT PREPARES THE INPUT DATA FOR THE SOLVER
% INPUTS OF FUNCTION 'INPUT_PERIODIC'
% ----------------------------------
% -------------------------------------------------------
%%%  2. Name of the file containing the mesh information
FUNinput.INPUTS.NameFileMesh = NameFileMesh ; %'mesh40k.msh' ;
FUNinput.INPUTS.NEUMANN_BOUNDARY_CONDITIONS = 'DISTRIBUTED_LOAD_Zmax' ; 'DISTRIBUTED_LOAD_Xmax_vert' ;  'DISTRIBUTED_LOAD_Xmax' ;   

% --------------------------------------------------------
% 3. MATERIAL PROPERTIES  (elastic)
% ----------------------
FUNinput.INPUTS.MATERIAL =  MATERIAL   ;
 
%----------------------------------------------
% ---- END INPUTS ----------------------------------

%----------------------------------------------------
% Calling Finite Element elastostatic program
% ---------------------------------------------------
%----------------------------------------------------
DATAOUT = FE_ELASTOSTATIC(FUNinput,DATA) ;
% ----------------------------------------------------


DATAOUT =   Strength_Laminate(DATAOUT,MATERIAL,DATA) ; 
