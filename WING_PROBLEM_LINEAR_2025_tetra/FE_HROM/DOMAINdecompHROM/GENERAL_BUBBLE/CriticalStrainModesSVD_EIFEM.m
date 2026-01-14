function [UpsilonDEF,S_upsilonDEF,V_upsilonDEF] =  CriticalStrainModesSVD_EIFEM(SNAPcompl,DATAoffline,MdomCHOL,PhiDEFbs)
%--------------------------------------------------------------------------
% [U, S_upsilonDEF, V_upsilonDEF] = CriticalStrainModesSVD_EIFEM(SNAPcompl, DATAoffline, MdomCHOL, PhiDEFbs)
%
% This function computes a basis of critical strain modes for nonlinear
% reduced-order modeling in the EIFEM framework. The function focuses on
% identifying modes associated with the most critical (i.e., unstable or
% highly energetic) deformation snapshots, distinguishing them from less
% relevant deformation states.
%
% It uses a **partitioned Singular Value Decomposition (SVD)** to separate
% basic, critical, and complementary modes, applying different tolerances to
% each block. The computation leverages the **Sequential Randomized SVD
% (SRSVD)** algorithm [Hernández et al., 2023] to efficiently handle blockwise
% matrices without requiring the full snapshot matrix to be built explicitly.
%
% INPUTS:
%   - SNAPcompl   : Cell array of complementary deformation snapshots from different subdomains or problems
%   - DATAoffline : Struct containing configuration data, especially:
%                    · IndexesCriticalSnapshots: cell array of vectors indicating critical snapshot indices
%                    · TOLSVD_complementary_modes: truncation tolerance for non-critical snapshots
%   - MdomCHOL    : Cholesky factor of the mass matrix Mdom (i.e., Mdom = MdomCHOL' * MdomCHOL)
%   - PhiDEFbs    : Matrix of basic deformation modes (e.g., affine or coarse-scale strain modes)
%
% OUTPUTS:
%   - U              : Mass-normalized basis of the computed modes
%   - S_upsilonDEF   : Diagonal matrix of singular values (mode importance)
%   - V_upsilonDEF   : Right singular vectors from the SVD decomposition
%
% FUNCTIONAL STEPS:
%   1. The function first checks consistency between the provided snapshot
%      data and the list of critical indices.
%
%   2. The matrix blocks for SVD are then assembled in the cell array A:
%       - A{1} contains the basic modes (normalized)
%       - For each entry in SNAPcompl:
%         - One block for critical snapshots (tolerance = 0 → fully retained)
%         - One block for complementary snapshots (tolerance > 0 → possibly truncated)
%
%   3. Each block is mass-normalized using MdomCHOL and rescaled with its
%      Frobenius norm to balance energy across partitions.
%
%   4. The partitioned SVD is computed using the **SRSVD algorithm**, which:
%       - Processes each block sequentially, computing the orthogonal
%         complement of the accumulated basis
%       - Applies block-specific truncation tolerances (RELTOL)
%       - Enables efficient memory usage and blockwise control of resolution
%         (see: J.A. Hernández et al., "CECM: A continuous empirical cubature
%         method with application to the dimensional hyperreduction of parameterized
%         finite element models", arXiv:2308.03877)
%
%   5. Finally, the mass-normalized basis U is mapped back to the physical
%      displacement space using the inverse of the Cholesky factor.
%
% NOTES:
%   - This function is especially useful in nonlinear problems with bifurcations,
%     instabilities, or path-dependent effects, where it is critical to capture
%     deformation mechanisms linked to specific loading events.
%   - The resulting modes are optimal (in least-squares sense) for representing
%     both critical and non-critical deformation behavior under a partitioned,
%     tolerance-controlled strategy.
% Joaquin A. Hernandez, 17-May-2025
%--------------------------------------------------------------------------

if nargin == 0
    load('tmp.mat')
end

 
if length(DATAoffline.IndexesCriticalSnapshots) ~= length(SNAPcompl)
    error('The list of critical snapshots should be provided for all tests (if empty, all snapshots are selected)')
end

% Partitioned SVD *********
% -------------------------
nproj = 1 + 2*length(SNAPcompl) ;
A  = cell(1,nproj) ;
RELTOL = zeros(nproj,1) ; 
% First block = basic modes
A{1} = MdomCHOL*PhiDEFbs ;
A{1} = A{1}/norm(A{1},'fro') ;  % Normalization
iproj_Acum = 1; 
TOL = DATAoffline.TOLSVD_complementary_modes ;

TOL_crit = DATAoffline.SVD_TOLERANCE_CriticalSnapshots ; 

for iprojLOC = 1:length(SNAPcompl)
    IndexesCriticalSnapshots_loc = DATAoffline.IndexesCriticalSnapshots{iprojLOC} ; 
    if  isempty(IndexesCriticalSnapshots_loc)
        IndexesCriticalSnapshots_loc = 1:size(SNAPcompl{iprojLOC},2) ; 
        IndCompl = [] ;
    else
        IndCompl = setdiff(1:size(SNAPcompl{iprojLOC},2),IndexesCriticalSnapshots_loc) ; 
    end 
    A{iproj_Acum+1}  =  MdomCHOL*SNAPcompl{iprojLOC}(:,IndexesCriticalSnapshots_loc) ;
     A{iproj_Acum+1}= A{iproj_Acum+1}/norm(A{iproj_Acum+1},'fro') ;
    RELTOL(iproj_Acum+1) = TOL_crit ; 
    iproj_Acum = iproj_Acum + 1; 
     A{iproj_Acum+1}  =  MdomCHOL*SNAPcompl{iprojLOC}(:,IndCompl) ;
    A{iproj_Acum+1}= A{iproj_Acum+1}/norm(A{iproj_Acum+1},'fro') ;
       RELTOL(iproj_Acum+1) = TOL ; 
    iproj_Acum = iproj_Acum + 1; 
end
 
DATAlocS = [] ;
[U,S_upsilonDEF,V_upsilonDEF] = SRSVD(A,RELTOL,DATAlocS) ;

UpsilonDEF = MdomCHOL\U ;  % All modes (including basic ones)