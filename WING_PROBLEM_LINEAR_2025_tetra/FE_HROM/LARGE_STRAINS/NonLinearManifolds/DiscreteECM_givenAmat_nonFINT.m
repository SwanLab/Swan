function [setPoints,wRED,ERROR_GLO,DATAOUT] = DiscreteECM_givenAmat_nonFINT(SNAPfint_lin,SNAPfint_non,DATA,wSTs,DATAoffline)
%--------------------------------------------------------------------------
% This is a modification of DiscreteECM_givenAmat, described below.
% We no longer apply blockwise the SVD to SNAPfint.
% August 25th 2025, Monday, Balmes 185, Barcelona.
%
% function [setPoints, wRED, ERROR_GLO, DATAOUT] = DiscreteECM_givenAmat(SNAPfint, DATA, wSTs, DATAoffline)
%
% Purpose:
%   Constructs a hyperreduced integration rule for nonlinear finite element
%   models using the Discrete Empirical Cubature Method (DECM).
%
%   The function receives a snapshot matrix of internal forces (SNAPfint),
%   performs orthogonalization to extract a reduced basis Q, and then applies
%   a greedy ECM algorithm to select a subset of integration points (Gauss points)
%   and their associated weights (wRED) that preserve the integrals over Q.
%
% Inputs:
%   - SNAPfint     : Snapshot matrix of internal forces (cell array or numeric)
%   - DATA         : Mesh and simulation data structure
%   - wSTs         : Vector of original Gauss point integration weights
%   - DATAoffline  : Structure with tolerance and algorithm settings
%       • .errorECM                      : target tolerance for cubature error
%       • .ListElementsExclude_fromGID  : (optional) file listing elements to exclude from candidate set
%       • .USE_SELECTIVE_DECM           : flag to choose ECM version
%
% Outputs:
%   - setPoints    : Indices of selected Gauss points for reduced integration
%   - wRED         : Reduced integration weights (corresponding to setPoints)
%   - ERROR_GLO    : Global integration error reported by the ECM routine
%   - DATAOUT      : Struct containing timing and diagnostics info
%
% Workflow:
%   1. Build orthogonal basis Q for the integrand from SNAPfint.
%   2. Exclude user-defined Gauss points (if specified in DATAoffline).
%   3. Apply DECM (Discrete Empirical Cubature Method) using Q and wSTs.
%   4. Compute and report actual relative integration error.
%
% Notes:
%   - Supports both original and updated ECM routines via USE_SELECTIVE_DECM flag.
%   - The integral error is validated post hoc by comparing exact and reduced integrals.
%   - If SNAPfint is a cell array, it's flattened for the final error computation.
%
% Author:
%   Joaquín A. Hernández, UPC, 12-JUL-2025
%  SEe  /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/112_NonLIN_ROM_RBF/02_1param_BEAM.mlx
%--------------------------------------------------------------------------
if nargin == 0
    load('tmp1.mat')
end


% Basis matrix for internal forces
% **********************************
%if iscell(SNAPfint)
SNAPfint_lin = cell2mat(SNAPfint_lin) ;
SNAPfint_non = cell2mat(SNAPfint_non) ;
%end
TIMEq = tic ;
[Q,SNAPfint] = QbasisMatrixIntegrand_givenSNAP_inelEL(SNAPfint_lin,SNAPfint_non,DATA,wSTs,DATAoffline) ;
DATAOUT.TIMEq = toc(TIMEq) ;


disp(['Number of internal force modes = ',num2str(size(Q,2))])


disp(['Time to assembly orthogonal matrix of the integrand =',num2str(DATAOUT.TIMEq)])



% Empirical cubature method
% -------------------------

%[setPoints,wRED,ERROR_GLO,DATAOUT]= EmpiricalCubatureMethod_orig(Q,[],wSTs,DATA_ECM)  ;

% CANDIDATE POINTS EXCLUDED FROM INFO IN FILE
% -------------------------------------------
TIMEdecm = tic ;

DATAoffline = DefaultField(DATAoffline,'ListElementsExclude_fromGID',[]) ;
if ~isempty(DATAoffline.ListElementsExclude_fromGID)
    ListElementsToExclude = load(DATAoffline.ListElementsExclude_fromGID) ;
    ListElementsToExclude = ListElementsToExclude(:,1) ;
    ngausELEM = DATA.MESH.ngaus_STRESS;
    ListGaussToExclude = small2large(ListElementsToExclude,ngausELEM) ;
    INDSEL  = setdiff(1:size(Q,1),ListGaussToExclude) ;
else
    INDSEL  = 1:size(Q,1) ;
end

DATA_ECM.IND_POINTS_CANDIDATES = INDSEL ;
DATAoffline = DefaultField(DATAoffline,'errorECM',0) ;
DATA_ECM.TOL = DATAoffline.errorECM ;
DATAoffline = DefaultField(DATAoffline,'USE_SELECTIVE_DECM',1) ; % = 1;
if DATAoffline.USE_SELECTIVE_DECM == 0
    % Version before 3-DEc-2021
    [setPoints,wRED,ERROR_GLO,DATAOUT]= EmpiricalCubatureMethod_orig(Q,[],wSTs,DATA_ECM)  ;
else
    % New version (after 3-DEc-2021)
    [setPoints,wRED,ERROR_GLO,DATAOUT]= DiscreteEmpiricalCubatureMethod(Q',wSTs,DATA_ECM)  ;
    
end

DATAOUT.TIMEdecm = toc(TIMEdecm) ;

disp(['Time to select the integration points =',num2str(DATAOUT.TIMEdecm)])


SNAPfint = cell2mat(SNAPfint) ;

IntegralExact =SNAPfint'*wSTs  ;
IntegrationError = IntegralExact - SNAPfint(setPoints,:)'*wRED ;
IntegrationError = norm(IntegrationError)/norm(IntegralExact);
disp(['Actual integration error using DECM= ',num2str(IntegrationError*100),' % (prescribed tolerance fint =',num2str(DATAoffline.errorFINT*100), '%'])