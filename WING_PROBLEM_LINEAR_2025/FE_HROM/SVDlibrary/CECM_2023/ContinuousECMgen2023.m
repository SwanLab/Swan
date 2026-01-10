function [ECMdata,DATA_OUT]= ContinuousECMgen2023(A,xFE,wFE,DATA,MESH,DATA_ECM)
%--------------------------------------------------------------------------
% function [ECMdata, DATA_OUT] = ContinuousECMgen2023(A,xFE,wFE,DATA,MESH,DATA_ECM)
%
% PURPOSE:
%   Implements the Continuous Empirical Cubature Method (CECM), which extends 
%   the Discrete Empirical Cubature Method (DECM) by introducing a sparsification
%   procedure for reducing the number of integration points while ensuring 
%   accurate integration of reduced basis functions. 
%
%   The algorithm constructs an orthogonal basis for the integrand functions 
%   (column space of A), computes a DECM rule, and optionally continues with 
%   the CECM step to yield a reduced quadrature rule with fewer integration points.
%
% INPUTS:
%   - A           : [N_gauss_points x N_samples] Sampling matrix (integrand values)
%   - xFE         : [N_gauss_points x d] Gauss point coordinates of the FE mesh
%   - wFE         : [N_gauss_points x 1] Gauss weights (including Jacobian)
%   - DATA        : Structure with general settings (e.g., exact integral)
%   - MESH        : Mesh structure containing connectivity and geometry
%   - DATA_ECM    : Structure with ECM/CECM options and parameters
%
% OUTPUTS:
%   - ECMdata     : Structure containing DECM/CECM quadrature data
%                   Fields include: 
%                      xCECM, wCECM (CECM points and weights),
%                      xDECM, wDECM (DECM points and weights),
%                      setElements (selected elements),
%                      Ninterpolation, INTexac (exact integral value)
%
%   - DATA_OUT    : Structure with extended outputs, including:
%                      - DATA_ECM    : Updated ECM control structure
%                      - MESH        : FE mesh
%                      - AUXVAR      : Auxiliary variables used in CECM
%                      - VAR_SMOOTH_FE: Interpolated integrand data
%                      - DATAOUTdecm : DECM history (selection and weights)
%                      - POLYINFO    : Info on polynomial interpolation elements
%
% THEORY:
%   The method uses a two-stage approach:
%   1. **Discrete ECM (DECM)**:
%      Constructs an interpolatory rule for a set of orthogonal basis functions
%      obtained via a weighted SVD of A. The DECM selects integration points among
%      the Gauss points and computes positive weights ensuring exact integration.
%
%   2. **Continuous ECM (CECM)**:
%      Starting from DECM points, a sparsification step removes points iteratively
%      by solving nonlinear constraints over polynomially interpolated integrand 
%      functions. The process aims to preserve integration accuracy while reducing
%      the number of points.
%
%   See Hernandez et al. (2024), *Computer Methods in Applied Mechanics and Engineering*,
%   vol. 418, article 116552.
%
% REFERENCE CODE:
%   - Based on scripts and ideas from:
%     https://github.com/Rbravo555/CECM-continuous-empirical-cubature-method
%   - Original test script: README_08.mlx
%
% DEPENDENCIES:
%   - GetBasisMatrixECM2: Computes weighted orthogonal basis for A
%   - DECMgeneral2023   : Performs greedy selection of DECM points
%   - ContinuousECM2023 : Executes sparsification for CECM
%   - DataPreparationCECM: Builds interpolants for continuous formulation
%   - VariousPLOTS_CECM2023: (Optional) Visualization of results
%   - PlotOrthogonalFunctToIntegrateCECM: Displays basis functions
%
% AUTHOR:
%   J.A. Hernández (UPC/CIMNE), adapted April–May 2024.
% Comments by ChatGPT4, 29-May-2024
%--------------------------------------------------------------------------



% Adaptation of ContinuousECM.m, general functions
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/ContinuousEmpiricalCubatureM/Paper_hernandez2021ecm/08_BACK_TO_POLY/README_08.mlx
if nargin == 0
    load('tmp3.mat')
end
ECMdata = [] ;
%DATA_ECM.NameFileMesh_ECM = [DATA.BASE_FOLDER,'DECMpoints'];
%DATA_ECM.NameFileMesh_CECM = [DATA.BASE_FOLDER,'CECMpoints'] ;
% Basis matrix for the column space of A
% ------------------------------------------------------------------
DATA_ECM = DefaultField(DATA_ECM,'Method_Evaluation_Basis_Integrand',1) ;  % % 1 --> FE fitting   / 2 ---> FE analytical evaluation
DATA_ECM = DefaultField(DATA_ECM,'NumberOfCECMpoints',[]) ;  % % 1 --> 
 if    isempty(DATA_ECM.NumberOfCECMpoints)
     DATA_ECM.NumberOfCECMpoints = 1e20 ;   % Apr-28-2024
 end

DATA = DefaultField(DATA,'ExactIntegral',[]) ;
if isempty(DATA.ExactIntegral)
    error('The exact integral should be calculated a priori')
end
[BasisA_squareW,VSinv] = GetBasisMatrixECM2(A,wFE,DATA_ECM) ;
%% DISCRETE EMPIRICAL CUBATURE METHOD (DECM)
% -------------------------------------------------------------------------
%[ECMdata,HYPERREDUCED_VARIABLES,DATAOUTdecm] = DECMgeneral(A,BasisA_squareW,wFE,xFE,DATA_ECM,MESH,DATA) ;
[ECMdata,HYPERREDUCED_VARIABLES,DATAOUTdecm] = DECMgeneral2023(A,BasisA_squareW,wFE,xFE,DATA_ECM,MESH,DATA) ;

% PLOT ORTHOGONAL FUNCTIONS TO BE INTEGRATED
DATA_ECM = DefaultField(DATA_ECM,'PLOT_ORTHOGONAL_FUNCTIONS_TO_BE_INTEGRATED',1) ;
PlotOrthogonalFunctToIntegrateCECM(HYPERREDUCED_VARIABLES.PHI,MESH,DATA_ECM) ; 

 


ECMdata.HISTORY.xDECM = DATAOUTdecm.HistoryPoints ;
ECMdata.HISTORY.wDECM = DATAOUTdecm.HistoryWeights ;
% ---------------------------------------------------------------------------
% DATA PREPARATION FOR Continuous EMPIRICAL CUBATURE METHOD (CECM)
% ---------------------------
DATA_ECM = DefaultField(DATA_ECM,'UseDECMpoints',0) ;
 
if DATA_ECM.UseDECMpoints == 0
    [VAR_SMOOTH_FE,DATA,DATA_ECM] = DataPreparationCECM(MESH,DATA,DATA_ECM,wFE,xFE,VSinv) ;
    
    % CONTINUOS EMPIRICAL CUBATURE METHOD (CECM)
    % -----------------------------------------------------------------------------
    [xNEW,wNEW,POLYINFO ,VAR_SMOOTH_FE,Ninterpolation,DATA_ECM,AUXVAR ] = ContinuousECM2023...
        (wFE,xFE,HYPERREDUCED_VARIABLES,DATA_ECM,VAR_SMOOTH_FE)  ;
    % --------------------------------------------------------------------------------------------
    
    %%%  PREPARING OUTPUT/PLOTTING FACILITIES
    % --------------------------------------------
    ECMdata.xCECM = xNEW ;
    ECMdata.wCECM = wNEW ;
    ECMdata.setElements = POLYINFO.setElements ; %  ELEMENTS_xNEW ;
    ECMdata.Ninterpolation = Ninterpolation ;
    ECMdata.INTexac = DATA.ExactIntegral;
    AUXVAR.DATAOUTdecm = DATAOUTdecm;
    
    VariousPLOTS_CECM2023(DATA_ECM,ECMdata,MESH,AUXVAR,VAR_SMOOTH_FE)  ;
    
    DATA_OUT.DATA_ECM = DATA_ECM ;
    DATA_OUT.ECMdata = ECMdata ;
    DATA_OUT.MESH = MESH ;
    DATA_OUT.AUXVAR = AUXVAR ;
    DATA_OUT.VAR_SMOOTH_FE = VAR_SMOOTH_FE ;
    DATA_OUT.DATAOUTdecm = DATAOUTdecm ;
    DATA_OUT.POLYINFO = POLYINFO ;
    
    DATA_OUT.indexPoints_DECM = HYPERREDUCED_VARIABLES.setPoints;
    
else
    
    ECMdata.xCECM =[];
    ECMdata.wCECM = [];
    
    if ~isempty(DATA_ECM.NumberOfCECMpoints) &&   length(HYPERREDUCED_VARIABLES.setElements) >  DATA_ECM.NumberOfCECMpoints     % JAHO, 21-Apr-2024
        disp('WARNING --> Number of integration points set by the user ...NumberOfCECMpoints')
        Npoints_selected = min(DATA_ECM.NumberOfCECMpoints,length(HYPERREDUCED_VARIABLES.setElements)) 
        if Npoints_selected ~= DATA_ECM.NumberOfCECMpoints
            disp(['Number of DECM points below the prescribed value '])
        end
     else
        Npoints_selected = length(HYPERREDUCED_VARIABLES.setElements) ; 
    end
    
    
    ECMdata.setElements = HYPERREDUCED_VARIABLES.setElements(1:Npoints_selected) ; %  ELEMENTS_xNEW ;
    ECMdata.Ninterpolation = [] ;
  %  ECMdata.INTexac = DATA.ExactIntegral;
  %  AUXVAR.DATAOUTdecm = DATAOUTdecm;
    
  %  VariousPLOTS_CECM2023(DATA_ECM,ECMdata,MESH,AUXVAR,VAR_SMOOTH_FE)  ;
    
    DATA_OUT.DATA_ECM = DATA_ECM ;
    ECMdata.wDECM = ECMdata.wDECM(1:Npoints_selected) ; % JAHO 21-Apr-2024
    ECMdata.xDECM = ECMdata.xDECM(1:Npoints_selected,:) ; 
    
    DATA_OUT.ECMdata = ECMdata ;
    DATA_OUT.MESH = MESH ;
    DATA_OUT.AUXVAR = [] ;
    DATA_OUT.VAR_SMOOTH_FE = [] ;
    DATA_OUT.DATAOUTdecm = [] ;
    DATA_OUT.POLYINFO = [] ;
    
    DATA_OUT.indexPoints_DECM = HYPERREDUCED_VARIABLES.setPoints(1:Npoints_selected);
    
end




