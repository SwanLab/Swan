function [stressMACRO DATAOUT DATA]= LaminateSTRESScal(strainINP,NameFileMesh,MATERIAL,DATA)
format long g
% Compute average generalized stresses (including moments) in a
% laminate  using  the extension to laminates of the
%   computational homogenization theory (by the finite element method)
%  Technical University of Catalonia, Nov-2015
%  % Joaquín A. HERNANDEZ, jhortega@cimne.upc.edu
% --------------------------------------------------------------------------------------------

if nargin ==0
    % INPUT DATA (default)
    % --- Generalized strains (input)
    strainINP = zeros(8,1) ;
    strainINP(1) = 0;   % ex
    strainINP(2) = 0 ;   % ey
    strainINP(3) = 0;   % exy
    strainINP(4) = 1/80 ;   %  Curvature k_x
    strainINP(5) = 0 ;   % Curvature k_y
    strainINP(6) = 0; %  Curvature k_xy
    strainINP(7) = 0 ; % shear gamma_yz
    strainINP(8) = 0 ; % shear gamma_xz
    % - Gid's file containing mesh data (COORDINATES, CONNECTIVITIES, MATERIAL and BOUNDARIES)
    NameFileMesh = 'lam10lay_8000.msh'; %'lam10lay_4RVE_7000.msh' ;'lam10lay_20000.msh';  'lam10lay_128.msh';'layer10_16rve.msh';  ; 'lam10lay_128.msh'  ;'lam10lay_4RVE_20000.msh';'lam10lay_14000.msh';'lam10lay_4RVE_7000.msh' ;'layer10_4rve.msh';   'lam10lay_8000.msh';'lam10lay_128.msh'  ;  'layer10_4rve_32000.msh'  ; ;'layer10_4rve_32000.msh' ; 'layer10_4rve_32000.msh' ; 'lam10lay_8000.msh' ; 'layer10_fineMESHthickness.msh';'lam10lay_8000.msh' ;
    NameWS_unitcell ='DATAWS/Celas_ejem1.msh.mat' ;  % 'DATAWS/Celas_iSO.msh.mat' ;
    NameWS_failureE = 'DATAWS/FAILURE_ejem1.msh.mat' ;
    % Elasticity matrix for each  layer (and angle of rotation)
    % The name of the .mat containg such information is to be provided  .
    % -----------------------
    %
    nplies = 8 ;
    ANG_FIB = [0 90 0 90 90 0 90 0] ;
    for iply = 1:nplies
        MATERIAL.PLY(iply).NAMEWS = NameWS_unitcell;  % Store the Celas in this folder
        MATERIAL.PLY(iply).ANGLE = ANG_FIB(iply) ;  % Angle subtended by the fiber and the x-axis
        MATERIAL.PLY(iply).NAMEWS_FAIL = NameWS_failureE ;
    end
    
    DATA.RECALCULATE_STIFFNESS =1;   % To avoid computing again the stiffness matrix (if =0) when computing
    % the average stresses for different input
    % strains, yet similar
    % material/geometric properties
    DATA.BOUNDARY_CONDlam = 'PERIODIC_uv_TOP_BOTTOM' ;    'PERIODIC'; 'ZERO_minTB' ;  'PERIODIC_minTOPB' ;  'PERIODIC_minTOPBzero' ;      'ZERO' ;   %PERIODIC';     % PERIODIC/ZERO
    DATA.ORDER_SHEAR_THEORY = 1 ;
    DATA.niterCONJG = 3000 ; % Number of iterations  conjugated gradient
    DATA.CALCULATE_averageSTRESS =2 ; %
    DATA.CALCULATE_STRENGTH_LAM = 1; 
    DATA.NUMBER_ELEMENTS_FAIL = 10;
end


% INPUT DATA
% ----------
FUNinput.NAME = 'INPUT_LAMINATES' ; % NAME OF THE FUNCTION THAT PREPARES THE INPUT DATA FOR THE SOLVER
% INPUTS OF FUNCTION 'INPUT_PERIODIC'
% ----------------------------------
% 1. MACROSCOPIC STRAINS
% -----------------------

FUNinput.INPUTS.strainINP =  strainINP ;
% -------------------------------------------------------
%%%  2. Name of the file containing the mesh information
FUNinput.INPUTS.NameFileMesh = NameFileMesh ; %'mesh40k.msh' ;
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
[DATAOUT DATA]= FE_ELASTOSTATIC(FUNinput,DATA) ;
% ----------------------------------------------------
% ----------------------------------------------------
% OUTPUT: AVERAGE STRESS ("MACROSCOPIC STRESS")
% -----------------------------------------------

stressMACRO = DATAOUT.stressAVG ;
disp('--------------------------------------------')
disp('Macroscopic stress')
disp(num2str(stressMACRO))



if isfield(DATA,'CALCULATE_STRENGTH_LAM') && DATA.CALCULATE_STRENGTH_LAM==1
    % ---------------------------------
    % Computing  ULTIMATE STRENGHT AND STRAINS4
    % ---------------------------------
    DATAOUT =   Strength_Laminate(DATAOUT,MATERIAL,DATA) ;
    %%%%
   
    
end
