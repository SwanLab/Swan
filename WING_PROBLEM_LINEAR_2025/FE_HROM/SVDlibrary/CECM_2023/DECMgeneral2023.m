function [ECMdata,HYPERREDUCED_VARIABLES,DATAOUTdecm] = DECMgeneral2023(A,U,wFE,xFE,DATA_ECM,MESH,DATA)
%--------------------------------------------------------------------------
% function [ECMdata, HYPERREDUCED_VARIABLES, DATAOUTdecm] = DECMgeneral2023(A, U, wFE, xFE, DATA_ECM, MESH, DATA)
%
% PURPOSE:
%   Implements the Discrete Empirical Cubature Method (DECM) for computing
%   an efficient quadrature rule with positive weights and reduced number
%   of Gauss integration points. The selected points and weights integrate
%   a given set of orthonormal basis functions (U) derived from a snapshot
%   matrix (A), and are used in reduced-order or hyperreduced FE models.
%
% INPUTS:
%   - A         : [M x N] Snapshot matrix of integrand functions (for error check)
%   - U         : [M x r] Orthonormal basis matrix (left singular vectors of A*sqrt(W))
%   - wFE       : [M x 1] Vector of Gauss weights (including Jacobian)
%   - xFE       : [M x d] Coordinates of Gauss points
%   - DATA_ECM  : Structure with DECM configuration options, including:
%                  * errorDECM (tolerance for integration error)
%                  * ListElementsToExclude (optional element exclusion)
%                  * IND_POINTS_CANDIDATES (candidate indices)
%   - MESH      : Structure containing mesh topology and ngausE (Gauss points per element)
%   - DATA      : General data structure; used here for optional error checking
%
% OUTPUTS:
%   - ECMdata   : Structure with DECM quadrature information:
%                  * xDECM, wDECM: coordinates and weights of selected integration points
%   - HYPERREDUCED_VARIABLES : Structure for use in hyperreduction:
%                  * setPoints     : indices of selected Gauss points
%                  * setElements   : corresponding FE elements (with repetition)
%                  * WdomRED       : DECM weights
%                  * PHI           : DECM basis matrix rescaled by 1/sqrt(wFE)
%   - DATAOUTdecm : Additional output from DiscreteEmpiricalCubatureMethod
%
% METHOD:
%   1. Preprocess candidate points: optionally exclude some elements from selection
%   2. Call `DiscreteEmpiricalCubatureMethod` to perform greedy selection of
%      quadrature points (among `INDSEL`) for the given basis U.
%   3. Map selected Gauss points to their parent elements (with/without repetition).
%   4. Assemble output variables and optionally compute integral errors (if requested).
%
% REMARKS:
%   - The DECM solves the linear system U^T W ≈ U^T W_REDUCED by selecting a minimal
%     set of rows (Gauss points) and associated positive weights.
%   - The resulting integration rule ensures accurate integration of the basis
%     functions (columns of U), with significant reduction in computational cost.
%   - Mapping from global Gauss points to elements is done using utility functions:
%       * small2large, large2smallREP, large2smallINCLUDEREP
%   - ErrorApproximateIntegral2 compares full and reduced integrals for A.
%
% REFERENCES:
%   - Hernández et al., "CECM: A continuous empirical cubature method with application
%     to the dimensional hyperreduction of parameterized finite element models,"
%     Comput. Methods Appl. Mech. Engrg., Vol. 418, 2024, 116552.
%
% SEE ALSO:
%   - DiscreteEmpiricalCubatureMethod
%   - GetBasisMatrixECM2
%   - ContinuousECMgen2023
%   - ErrorApproximateIntegral2
%
% AUTHOR:
%   J.A. Hernández (UPC/CIMNE), DECM adaptation for 2023 hyperreduction framework.
% Comments by ChatGPT4, 29-May-2025
%--------------------------------------------------------------------------

if nargin  == 0
    load('tmp1.mat')
end

% List of elements to be excluded from the initial SET
DATA_ECM = DefaultField(DATA_ECM,'ListElementsToExclude',[]) ;

INDSEL = 1:length(wFE) ;

DATA_ECM = DefaultField(DATA_ECM,'IND_POINTS_CANDIDATES',INDSEL) ; 
 

if ~isempty(DATA_ECM.ListElementsToExclude)
    ListGaussToExclude = small2large(DATA_ECM.ListElementsToExclude,MESH.ngausE) ;
    INDSEL  = setdiff(DATA_ECM.IND_POINTS_CANDIDATES,ListGaussToExclude) ;
else
    INDSEL = DATA_ECM.IND_POINTS_CANDIDATES ; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discrete Empirical cubature method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DATA_ECM = DefaultField(DATA_ECM,'errorDECM',0) ;
DATA_DECM = [] ;
DATA_DECM.TOL = DATA_ECM.errorDECM ;
DATA_DECM.IND_POINTS_CANDIDATES = INDSEL ;
DATA_DECM.STORE_INFO_ITERATIONS = 1; 
[zDECM,wDECM,~,DATAOUTdecm]= DiscreteEmpiricalCubatureMethod(U',wFE,DATA_DECM)  ;
ECMdata.xDECM = xFE(zDECM,:);
ECMdata.wDECM = wDECM ;
% Determining the indices of the associated elements
setElements = large2smallREP(zDECM,MESH.ngausE) ;
setElementsREP = large2smallINCLUDEREP(zDECM,MESH.ngausE) ;

disp('****************************+')
disp(['List of selected m = ',num2str(length(setElements)),' elements'])
disp(num2str(setElements'))
%clipboard('copy',num2str(setElements'));

HYPERREDUCED_VARIABLES.setPoints = zDECM ;  % SEt integration points
HYPERREDUCED_VARIABLES.setElements = setElementsREP ;  % Set associated elements
HYPERREDUCED_VARIABLES.WdomRED = wDECM ;  % Set associated WEights
HYPERREDUCED_VARIABLES.PHI = bsxfun(@times,U,1./sqrt(wFE)) ;


% ------------------------------------------------
% Computing integration errors due to the DECM
% -----------------------------------------------

DATA = DefaultField(DATA,'CalculateErrorIntegralSnapshotMatrices',1) ; % = 0 ;

if DATA.CalculateErrorIntegralSnapshotMatrices == 1
        ErrorApproximateIntegral2(A,HYPERREDUCED_VARIABLES.PHI,wFE,DATA_ECM,zDECM,wDECM,DATA)  ;
end
