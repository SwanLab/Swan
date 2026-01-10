function Train1dom_all(DATAlauch,INPUT_PROBLEMS,NameWS_EIFE_oper,DATAoffline)
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
% 30-JAN-2023/25-Feb-2023/17-Nov-2023, JAHO
% ------------------------------------------------------
% Inputs
% ----------------
if nargin == 0
   load('tmp.mat')     
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
    OFFLINE_STEPSfun(INPUT_PROBLEMS,DATAoffline) ;    
end
% EMPIRICAL INTERSCALE FINITE ELEMENT OPERATORS (EIFE)
DATAlauch = DefaultField(DATAlauch,'STORE_additional_FE_information',0) ; 
if DATAlauch.ONLINE_OPERATORS == 1
    [EIFEoper,OPERFE,MESHdom,GEOMATRICES,MODES] = EIFE_operatorsFUN(INPUT_PROBLEMS,DATAoffline)  ;
    if DATAlauch.STORE_additional_FE_information == 1
        save(NameWS_EIFE_oper,'EIFEoper','OPERFE','MESHdom','GEOMATRICES','MODES') ;
    else
        save(NameWS_EIFE_oper,'EIFEoper') ; 
    end
end
