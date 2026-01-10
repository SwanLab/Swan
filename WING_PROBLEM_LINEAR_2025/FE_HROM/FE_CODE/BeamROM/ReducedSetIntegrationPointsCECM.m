function [HYPERREDUCED_VARIABLES,MSG ] =...
    ReducedSetIntegrationPointsCECM(DATAIN,DATAROM,BASES,DATA_REFMESH,MSG)

if nargin == 0
    load('CheckCECM.mat')
    
end
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/ContinuousEmpiricalCubatureM/Paper_hernandez2021ecm/09_ASSESSMENTpaper/MultiscaleHROM/READMEmulti.mlx

DATAIN = DefaultField(DATAIN,'ISNONLINEAR',0) ;
DATAIN = DefaultField(DATAIN,'CUBATURE',[]) ;

 % Continuous
[HYPERREDUCED_VARIABLES,MSGloc,Wdom] = CECM_multi(DATAIN,DATAROM,BASES,DATA_REFMESH) ;

% 
% %HYPERREDUCED_VARIABLES.setPoints = setPoints ;  % SEt integration points
% HYPERREDUCED_VARIABLES.setElements = setElements ;  % Set associated elements
HYPERREDUCED_VARIABLES.WdomRED = HYPERREDUCED_VARIABLES.wCECM ;  % Set associated WEights
% %HYPERREDUCED_VARIABLES.Celas_Bdom = Celas_Bdom ;  % Product Elastic Matrix times  B matrix
% HYPERREDUCED_VARIABLES.nstrain = nstrain ;
%HYPERREDUCED_VARIABLES.CelasBdomALL = CelasBdomALL ;

% Discrete
%[HYPERREDUCED_VARIABLES,MSGloc,Wdom] = EmpiricalCubatureMethod_points(DATAIN,DATAROM,BASES,DATA_REFMESH) ;


for iii = 1:length(MSGloc)
    MSG{end+1} = MSGloc{iii} ;
end

% Reconstruction matrix
% ---------------------
DATAIN = DefaultField(DATAIN,'COMPUTE_RECONSTRUCTION_MATRIX',0) ;
if  DATAIN.COMPUTE_RECONSTRUCTION_MATRIX  == 1 || DATAIN.ISNONLINEAR==1
    error('Option not implemented yet')
    %     BasisS = BASES.STRESSES.U  ;
    %     setGauss = small2large(HYPERREDUCED_VARIABLES.setPoints,HYPERREDUCED_VARIABLES.nstrain) ;
    %     coeff = (BasisS(setGauss,:)'*BasisS(setGauss,:))\BasisS(setGauss,:)' ;
    %     HYPERREDUCED_VARIABLES.ReconsStresses = BasisS*coeff ;
    %
    %     % HOMOGENIZATION OPERATOR  (19-JAN-2020, for the paper)
    %     % --------------------------------------------------------------------------------------------
    %     load(DATAIN.NAME_WS_MODES,'Wdom')
    %     nstrain = size(BasisS,1)/length(Wdom);
    %     HOMOGENIZATION_OPERATOR = zeros(nstrain,size(HYPERREDUCED_VARIABLES.ReconsStresses,2)) ;
    %
    %     for istrain = 1:nstrain
    %         COMP = istrain:nstrain:size(HYPERREDUCED_VARIABLES.ReconsStresses,1) ;
    %         HOMOGENIZATION_OPERATOR(istrain,:) = Wdom'*HYPERREDUCED_VARIABLES.ReconsStresses(COMP,:) ;
    %     end
    %     HYPERREDUCED_VARIABLES.HOMOGENIZATION_OPERATOR = HOMOGENIZATION_OPERATOR ;
    %     % -------------------------------------------------------------------------------------------
    
end

%--------------------------------------
MSG{end+1} = '-----------------------------------------' ;
MSG{end+1} = 'Position ECM points (ans weights )' ;
MSG{end+1} = '-----------------------------------------' ;

DATAIN = PrintingGid_ECMpoints(DATAIN,DATA_REFMESH,HYPERREDUCED_VARIABLES) ;

for iprint =  1:length(DATAIN.MSGPRINT)
    MSG{end+1} = DATAIN.MSGPRINT{iprint} ;
end

% -------------------------------------
MSG{end+1} = '-----------------------------------------' ;



end

