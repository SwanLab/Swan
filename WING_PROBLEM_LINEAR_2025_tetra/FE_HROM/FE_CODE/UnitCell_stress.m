function stressMACRO = UnitCell_stress(strainMACRO,NameFileMesh,MATERIAL,RECALCULATE_STIFFNESS)
format long g
% Compute average stresses in a carbon-fiber/matrix unit cell
%  Technical University of Catalonia, Nov-2015
%  % Joaquín A. HERNANDEZ, jhortega@cimne.upc.edu
% --------------------------------------------------------------------------------------------

if nargin ==0
    % INPUT DATA (default)
    % --------------------------------
    strainMACRO = zeros(6,1) ;
    strainMACRO(1) = 1;   % ex
    strainMACRO(2) = 0 ;   % ey
    strainMACRO(3) = 0;   % ez
    strainMACRO(4) = 0 ;   % gamma_yz
    strainMACRO(5) = 0 ;   % gamma_xz
    strainMACRO(6) = 0; % gamma_xy
    
    NameFileMesh = 'mesh40k.msh' ;
    NameFileMesh = 'mesh_090vf.msh' ;
    % Transverse isotropic material
    MATERIAL.FIBER.E(1) = 230e3 ; %MPa  % Longitudinal elastic modulus
    MATERIAL.FIBER.E(2) = 8e3 ; %MPa   % Transverse elastic modulus
    MATERIAL.FIBER.nu(1) = 0.25 ;  % Major poisson's ratio (contraction in the transverse directionupon an extension in the fiber direction)
    MATERIAL.FIBER.nu(2) = 0.3 ;  % Transverse poisson's ratio 
    MATERIAL.FIBER.Gshear(1) = 27.3e3 ; %MPa % In-plane shear modulus
    MATERIAL.FIBER.Gshear(2)= 3.08e3 ; % MPa % Transverse shear modulus
    MATERIAL.FIBER.INDEX = 1;  % Material index 
    %%%
    E = 3e3 ; nu = 0.3 ; 
    MATERIAL.MATRIX.E(1) = E  ; % MPa Elastic modulus   % MATRIX
    MATERIAL.MATRIX.nu(1) = nu ; % Poisson's ratio
    MATERIAL.MATRIX.Gshear(1) = E/2/(1+nu)  ;  % In-plane shear modulus
    MATERIAL.MATRIX.INDEX = 2; % Material index
    
    
    RECALCULATE_STIFFNESS = 1;
    DATA.REFERENCE_POINT = 'CENTER';  % Options: 'CENTER','CORNER' (FOR COMPUTING TOTAL DISPLACEMENTS )
    DATA.STORE_STIFFNESS = 0 ; %Store stiffness matrix 
    DATA.niterCONJG = 10000 ; % Number of iterations  conjugated gradient
    DATA.tolCONJG = 1e-10 ; % Tolerance
end


% INPUT DATA
% ----------
FUNinput.NAME = 'INPUT_PERIODIC' ; % NAME OF THE FUNCTION THAT PREPARES THE INPUT DATA FOR THE SOLVER
% INPUTS OF FUNCTION 'INPUT_PERIODIC'
% ----------------------------------
% 1. MACROSCOPIC STRAINS
% -----------------------

FUNinput.INPUTS.strainMACRO =  strainMACRO ;
% -------------------------------------------------------
%%%  2. Name of the file containing the mesh information
FUNinput.INPUTS.NameFileMesh = NameFileMesh ; %'mesh40k.msh' ;
% --------------------------------------------------------
% 3. MATERIAL PROPERTIES  (elastic)
% ----------------------


FUNinput.INPUTS.MATERIAL =  MATERIAL   ;
%----------------------------------------------
DATA.RECALCULATE_STIFFNESS = RECALCULATE_STIFFNESS ; % To avoid computing again the stiffness matrix (if =0) when computing
% the average stresses for different input
% strains, yet similar
% material/geometric properties
% ---- END INPUTS ----------------------------------

%----------------------------------------------------
% Calling Finite Element elastostatic program
% ---------------------------------------------------
%----------------------------------------------------
DATAOUT = FE_ELASTOSTATIC(FUNinput,DATA) ;
% ----------------------------------------------------
% ----------------------------------------------------
% OUTPUT: AVERAGE STRESS ("MACROSCOPIC STRESS")
% -----------------------------------------------

stressMACRO = DATAOUT.stressAVG ;
disp('--------------------------------------------')
disp('Macroscopic stress')
disp(num2str(stressMACRO))
