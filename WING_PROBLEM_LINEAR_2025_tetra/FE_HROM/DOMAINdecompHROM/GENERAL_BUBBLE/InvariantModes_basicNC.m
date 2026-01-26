function [PhiDEFbs_INV,MdomCHOL,DATAoffline] = InvariantModes_basicNC(SNAPbasic,DATAcommon,PhiRB,Mdom,DATAoffline,MESH,DATA,Kstiff,Mintf)
%--------------------------------------------------------------------------
% FUNCTION: InvariantModes_basicNC
%
% PURPOSE:
%   Extracts *invariant* (elastic) deformation modes from a set of basic
%   training displacements in the small strain regime, by filtering out
%   rigid body motions and performing a weighted SVD.
%   These invariant modes serve as the foundation for constructing reduced
%   order bases in multiscale methods (e.g., EIFEM).
%
% DESCRIPTION:
%   - Filters zero columns from input snapshots (due to inactive DOFs)
%   - Projects displacement snapshots orthogonally to the space of rigid body modes
%   - Applies a weighted SVD (optionally with Cholesky preconditioning)
%   - Selects invariant modes based on the singular value decay or a user-defined number
%   - Optionally, uses invariant modes provided by the user via another mesh
%   - Provides geometric checks to ensure mesh compatibility in such case
%   - Plots the resulting invariant modes and evaluates their conjugacy with interface reactions
%
% INPUTS:
%   - SNAPbasic   : Matrix of displacement snapshots (basic/elastic tests)
%   - DATAcommon  : Configuration structure (may include user-defined invariant modes)
%   - PhiRB       : Matrix of rigid body modes
%   - Mdom        : Domain geometric mass matrix
%   - DATAoffline : Offline processing settings (tolerances, switches)
%   - MESH        : Mesh of the current subdomain
%   - DATA        : Simulation configuration (includes geometry)
%   - Kstiff      : Global stiffness matrix of the subdomain
%   - Mintf       : Interface geometric mass matrix (for plotting)
%
% OUTPUTS:
%   - PhiDEFbs_INV: Matrix containing the selected invariant deformational modes
%   - MdomCHOL    : Cholesky factor of Mdom (used for subsequent orthogonalization)
%   - DATAoffline : Possibly updated with mode counts or configuration flags
%
% PROCEDURE:
%   1. Remove zero (inactive) DOFs from input.
%   2. Orthogonalize displacements w.r.t. rigid body modes.
%   3. Perform weighted SVD on the resulting displacements.
%   4. Analyze the decay of singular values to select invariant modes.
%   5. If specified, use user-provided invariant modes (with coordinate matching).
%   6. Compute singular values of both displacements and corresponding interface forces.
%   7. Plot results, and optionally plot contributions along user-defined edges.
%
% KEY OPTIONS (DATAoffline fields):
%   - TOLSVD_INVARIANT_BASIC_MODES           : Relative SVD tolerance for selecting modes
%   - NumberOfInvariantModes_SVD             : Max number of invariant modes to keep
%   - USE_CHOLESKY_DECOMPOSITION             : If true, use Cholesky preconditioning
%   - THRESHOLD_RATIO_CONSECUTIVE_SINGLE_VALUES_INVARIANT_DEF_MODES: Ratio threshold
%   - PLOT_INFO_EDGES                        : If true, highlights deformation over element edges
%
% SEE ALSO:
%   - SprojDEF_operator
%   - WSVDT, SVDT
%   - PlotModesDEF_SubdomainLevel
%   - PlotInvModesEDGES
%
% AUTHOR:
%   Joaquín A. Hernández Ortega (UPC/CIMNE)
%   Created: 10-Apr-2024, Barcelona
%   Comments by ChatGPT3, 13-May-2025
%--------------------------------------------------------------------------


% Determination of invariant modes, for elastic deformational modes
% JAHO, 10-Apr-2024, UPC, Barcelona
% Modification InvariantModes_basic.m, case no complementary modes
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/107_EIFEM_plast3D/02_BEAMQ8_IprofE.mlx
% -------------------------------
if nargin == 0
    load('tmp1.mat')
end

% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/107_EIFEM_plast3D/07_2020_method.mlx
DATAcommon = DefaultField(DATAcommon,'PhiDEFbs_GIVEN_OTHER_MESH',[])  ;
b = MESH.faceDOFSall;


% REMOVING ZERO COLUMNS, FOR ELASTIC TESTS
%
SNAPbasic_NORM = sum(SNAPbasic.^2,1) ;
nonzerosCOL = find(SNAPbasic_NORM~=0);
SNAPbasic = SNAPbasic(:,nonzerosCOL) ;

 


%%% STEP 2.  INVARIANT MODES
% -------------------------------------------------------------------
DATAcommon = DefaultField(DATAcommon,'INFO_BASIC_DEFORMATIONAL_MODES',[]) ;
DATAoffline = DefaultField(DATAoffline,'USE_CHOLESKY_DECOMPOSITION',1) ; % JAHO 22-Apr-2024
% To reduce offline cost
USE_CHOLESKY_DECOMPOSITION = DATAoffline.USE_CHOLESKY_DECOMPOSITION ;

% Purging rigid-body displacements (small strains, because the basic training tests are in  SMALL strains)
SNAPinv = SprojDEF_operator(PhiRB,Mdom,SNAPbasic) ;
% ----


%DATAoffline.TOLSVD_INVARIANT_BASIC_MODES = 1e-4 ;   % Determine invariant modes
DATAoffline = DefaultField(DATAoffline,'TOLSVD_INVARIANT_BASIC_MODES',1e-4) ;
DATAoffline = DefaultField(DATAoffline,'NumberOfInvariantModes_SVD',1e20) ;   



DATAlocc.TOL = DATAoffline.TOLSVD_INVARIANT_BASIC_MODES;
if  USE_CHOLESKY_DECOMPOSITION == 1
    [ PhiDEFbs,Sbs,Vbs,MdomCHOL] = WSVDT( SNAPinv,Mdom,DATAlocc) ;
else
    [ PhiDEFbs,Sbs,Vbs,MdomCHOL] = WSVDT( SNAPinv,[],DATAlocc) ;
end

% Filtering suspicious modes/modes given from other tests
 PhiDEFbs_INV = CheckSuspiciousDefMODES_EIFEM(DATAcommon,Sbs,DATAoffline,PhiDEFbs,MESH) ; 


DATA.NAME_BASE = 'Deformational_ALL';
PlotModesDEF_SubdomainLevel(DATA,PhiDEFbs_INV,MESH);
disp(['*****************************************************************************************'])


% COMPUTING ANGLES FORMED BY THE DEFORMATIONAL MODES AT THE BOUNDARY AND
% THE CORRESPONDING REACTIONS 
% -------------------------------------------------------------------
[UUUdef,SSS,VVV] = SVDT(PhiDEFbs_INV(b,:)) ;
disp(['Singular values   PhiDEFbs_INV(b,:)'])
SSS/SSS(1)
%
PsiSE_INV = Kstiff(b,:)*PhiDEFbs_INV ;
%
[PsiSE_INV,SSSse,VVVse] = SVDT(PsiSE_INV) ;
disp(['Cosine angles formed by invariant deformational modes and their conjugate interface forces'])
[Uconj,Sconj,Vconj ] = SVDT(PsiSE_INV'*UUUdef) ;
Sconj
%
% disp(['Deformational modes sorted according to the contribution to the work done by the reactive modes'])
%

% 20-April-2025
% See  /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/107_EIFEM_plast3D/09_BACK_to_HEXAG.mlx
DATA = DefaultField(DATA,'INFO_EDGES',[]); 
DATAoffline = DefaultField(DATAoffline,'PLOT_INFO_EDGES',0) ; % Up to May 5th 2025, only valid for Q9 elements
 if ~isempty(DATA.INFO_EDGES) && DATAoffline.PLOT_INFO_EDGES == 1
    PlotInvModesEDGES(DATA,PhiDEFbs_INV(b,:),Mintf)
end


