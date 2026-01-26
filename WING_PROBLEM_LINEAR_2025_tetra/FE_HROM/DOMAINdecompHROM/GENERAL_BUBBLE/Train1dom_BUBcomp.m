function Train1dom_BUBcomp(DATAlauch,INPUT_PROBLEMS,NameWS_EIFE_oper,DATAoffline)

%--------------------------------------------------------------------------
% FUNCTION: Train1dom_BUBcomp
%
% PURPOSE:
%   Perform the training and offline stages of the Empirical Interscale 
%   Finite Element Method (EIFEM) for a single domain with bubble-compressed
%   multiscale enrichment. Depending on user flags, the script can:
%
%     - Run full-field simulations for training (snapshots)
%     - Compute reduced-order basis functions and select empirical integration points
%     - Generate online reconstruction and integration operators
%
%   These steps produce a compressed set of operators encapsulated in the
%   structure `EIFEoper`, suitable for fast online evaluation of internal forces.
%
% INPUTS:
%   - DATAlauch: Structure with execution flags (TRAINING, OFFLINE_modes_points, ONLINE_OPERATORS)
%   - INPUT_PROBLEMS: Name of the function that defines the geometry, mesh,
%                     material parameters, and boundary conditions
%   - NameWS_EIFE_oper: Path where the output binary file will be stored
%   - DATAoffline: Structure defining the target subdomain and boundary (LABEL_DOMAIN, LABELS_FACES)
%
% OUTPUT:
%   - EIFEoper: Structure containing operators for online ROM evaluation, including:
%       * INTforces:   Reduced internal force operator
%       * BodyForces:  Reduced body force operator
%       * RECONSTRUCTION: Basis matrices and enrichment
%       * MESH:        Geometry and connectivity of reduced model
%
% LOCATION OF MAIN ROUTINES CALLED:
%   - Training_MULTI:     Full-field simulations to collect snapshots
%   - OFFLINE_STEPSbubCOMPL: SVD compression and empirical cubature
%   - EIFE_operatorsBUBgen: Generation of reduced operators (robust version)
%
% DEPENDENCIES:
%   - Requires MATLAB_CODES_FEHROM in the path
%
% AUTHOR:
%   Joaquín A. Hernández Ortega (UPC-CIMNE)
%   Created: 30-Jan-2023 — Updated: 25-Feb-2023
%   Comments by ChatGPT4, 13-May-2025
%--------------------------------------------------------------------------


% SEe /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/05_Assessment.mlx
% SCRIPT FOR RUNNING BOTH THE TRAINING STAGE AND THE OFFLINE STAGE OF A
% MULTISCALE ROM
% Inputs: INPUT_PROBLEMS: file definining the training domain, boundary conditions, mesh,
% and so on
% LABEL_DOMAIN =  Label of the domain we wish to study
% LABELS_FACES =  Label surfaces/lines definining the boundary of the
% domain
% OUTPUT: Data structure EIFEoper  (Empirical Interescale Finite Element
% Method) operators, stored in
% the binary file NameWS_EIFE_oper
% EIFEoper =
%
%   struct with fields:
%
%          INTforces: [1×1 struct]
%         BodyForces: [1×1 struct]
%     RECONSTRUCTION: [1×1 struct]
%               MESH: [1×1 struct]

% -----------------------------------------------------
% 30-JAN-2023/25-Feb-2023, JAHO
% ------------------------------------------------------
% Inputs
% ----------------
if nargin == 0
    INPUT_PROBLEMS ='DefineInputs_HOMOGq2' ; 'DefineInputs_HOMOGq' ;'DefineInputs_HOMOG' ; 'DefineInputs_HOMOGbig' ;   % 'DefineInputs_3DLIN'; 'DefineInputs_HOMOG_quad' ;
    % Offline
    DATAoffline.LABEL_DOMAIN = 1; % LABEL IDENTIFYING THE SURFACE/VOLUME WE WISH TO STUDY
    DATAoffline.LABELS_FACES = [7,8,5,6] ;  [6,7,8,5] ;   % LABELS IDENTIFYING THE BOUNDARY OF SUCH A DOMAIN
    NameWS_EIFE_oper  = [cd,filesep,'EIFE_LIBRARY',filesep,'QUADlinHOMOG_p03_2.mat'] ;  % Binary file where to store the  resulting
    % Empirical Interscale FE obthect
    
    %
    DATAlauch.TRAINING =0;
    DATAlauch.OFFLINE_modes_points =0;
    DATAlauch.ONLINE_OPERATORS= 1;
    
    
    DATAoffline.errorFINT  = 0 ;
    DATAoffline.errorMASSmatrixDECM  =0 ;
    DATAoffline.UseDECMpoints = 0;

    
end





if ~exist('Quadratic3N')
    addpath(genpath(getenv('MATLAB_CODES_FEHROM')))  ;
end


 
if DATAlauch.TRAINING ==1
    % TRAINING FE SIMULATIONS
    Training_MULTI(INPUT_PROBLEMS,DATAoffline) ;
end
if DATAlauch.OFFLINE_modes_points ==1
    % BASIS MATRICES, AND INTEGRATION POINTS VIA THE CECM
    OFFLINE_STEPSbubCOMPL(INPUT_PROBLEMS,DATAoffline) ;
     
    
end
% EMPIRICAL INTERSCALE FINITE ELEMENT OPERATORS (EIFE)

[AAA,BBB,CCC] =fileparts(NameWS_EIFE_oper) ; 
if ~exist(AAA)
    mkdir(AAA) ; 
end


if DATAlauch.ONLINE_OPERATORS == 1
    
    % 
    DATAlauch = DefaultField(DATAlauch,'EnableResidualBasedMethodForTestingPurposes',0) ; 
    if  DATAlauch.EnableResidualBasedMethodForTestingPurposes == 1 
        disp('WArning: Option abandoned on Feb 23th in favor of EIFE_operatorsBUBgen.m')
        disp('You should only enable this if you want to recover previous versions....')
    EIFEoper = EIFE_operatorsBUB(INPUT_PROBLEMS,DATAoffline)  ;
    else 
         EIFEoper = EIFE_operatorsBUBgen(INPUT_PROBLEMS,DATAoffline)  ;
    end
    
    
    save(NameWS_EIFE_oper,'EIFEoper') ;
end
