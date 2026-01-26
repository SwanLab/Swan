function [ECMdata,DATA_OUT]= ContinuousECMgen(A,xFE,wFE,DATA,MESH,DATA_ECM)
% Adaptation of ContinuousECM.m, general functions
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/ContinuousEmpiricalCubatureM/Paper_hernandez2021ecm/08_BACK_TO_POLY/README_08.mlx
if nargin == 0
    load('tmp1.mat')
end
ECMdata = [] ;
DATA_ECM.NameFileMesh_ECM = [DATA.BASE_FOLDER,'DECMpoints'];
DATA_ECM.NameFileMesh_CECM = [DATA.BASE_FOLDER,'CECMpoints'] ;
% Basis matrix for the column space of A
% ------------------------------------------------------------------
DATA_ECM = DefaultField(DATA_ECM,'Method_Evaluation_Basis_Integrand',1) ;  % % 1 --> FE fitting   / 2 ---> FE analytical evaluation


DATA = DefaultField(DATA,'ExactIntegral',[]) ;
if isempty(DATA.ExactIntegral)
    error('Adapt the code so that it can use GetBasisMatrixECM2 (more effective)')
    [U,INTexac,VSinv,A] = GetBasisMatrixECM(A,wFE,DATA_ECM) ;   % Before 25-Oct-2022
    % A is pre-multiplied by sqrt(W)
    DATA.ExactIntegral = INTexac;
else
    [U,VSinv] = GetBasisMatrixECM2(A,wFE,DATA_ECM) ;
    % INTexac = DATA.ExactIntegral ;
end

%% DISCRETE EMPIRICAL METHOD (and error computations, if any)
% -------------------------------------------------------------------------
[ECMdata,HYPERREDUCED_VARIABLES,DATAOUTdecm] = DECMgeneral(A,U,wFE,xFE,DATA_ECM,MESH,DATA) ;
% ---------------------------------------------------------------------------

% DATA PREPARATION FOR CECM
% ---------------------------
[VAR_SMOOTH_FE,DATA,DATA_ECM] = DataPreparationCECM(MESH,DATA,DATA_ECM,wFE,xFE,VSinv) ; 

% CONTINUOS EMPIRICAL CUBATURE METHOD (CECM)
% -----------------------------------------------------------------------------

[xNEW,wNEW,POLYINFO ,VAR_SMOOTH_FE,Ninterpolation,DATA_ECM,AUXVAR ] = ConstrainedCECM...
    (wFE,xFE,HYPERREDUCED_VARIABLES,DATA_ECM,VAR_SMOOTH_FE)  ;

% --------------------------------------------------------------------------------------------

%%%  PREPARING OUTPUT/PLOTTING FACILITIES
% --------------------------------------------
ECMdata.xCECM = xNEW ;
ECMdata.wCECM = wNEW ;
ECMdata.setElements = POLYINFO.ELEMENTS_CONTAINING_xNEW ; %  ELEMENTS_xNEW ;
ECMdata.Ninterpolation = Ninterpolation ;
ECMdata.INTexac = DATA.ExactIntegral;
AUXVAR.DATAOUTdecm = DATAOUTdecm; 

VariousPLOTS_CECM(DATA_ECM,ECMdata,MESH,AUXVAR,VAR_SMOOTH_FE)  ;

DATA_OUT.DATA_ECM = DATA_ECM ;
DATA_OUT.ECMdata = ECMdata ;
DATA_OUT.MESH = MESH ;
DATA_OUT.AUXVAR = AUXVAR ;
DATA_OUT.VAR_SMOOTH_FE = VAR_SMOOTH_FE ;
 DATA_OUT.DATAOUTdecm = DATAOUTdecm ;

DATA_OUT.indexPoints_DECM = HYPERREDUCED_VARIABLES.setPoints;  





