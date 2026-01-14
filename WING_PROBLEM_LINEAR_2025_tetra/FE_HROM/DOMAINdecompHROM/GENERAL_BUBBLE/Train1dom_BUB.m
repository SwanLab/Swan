function Train1dom_BUB(DATAlauch,INPUT_PROBLEMS,NameWS_EIFE_oper,DATAoffline)
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
    OFFLINE_STEPSbub(INPUT_PROBLEMS,DATAoffline) ;
     
    
end
% EMPIRICAL INTERSCALE FINITE ELEMENT OPERATORS (EIFE)
if DATAlauch.ONLINE_OPERATORS == 1
    EIFEoper = EIFE_operatorsBUB(INPUT_PROBLEMS,DATAoffline)  ;
    save(NameWS_EIFE_oper,'EIFEoper') ;
end
