function [UpsilonDEF,S_upsilonDEF,V_upsilonDEF] = StrainModesNonlinearSVD_EIFEM(SNAPcompl,DATAoffline,Mdom,PhiDEFbs,PhiRB,MdomCHOL,MESH)
%--------------------------------------------------------------------------
% [UpsilonDEF, S_upsilonDEF, V_upsilonDEF] = StrainModesNonlinearSVD_EIFEM(SNAPcompl, DATAoffline, Mdom, PhiDEFbs, PhiRB, MdomCHOL)
%
% This function computes a set of strain-related deformation modes using a
% partitioned Singular Value Decomposition (SVD) tailored for nonlinear
% reduced-order modeling within the EIFEM (Empirical Interscale Finite
% Element Method) framework. The input snapshots are orthogonalized against
% rigid body motions and mass-normalized prior to applying SVD.
%
% INPUTS:
%   - SNAPcompl   : Cell array of complementary deformation snapshots
%   - DATAoffline : Structure containing offline configuration parameters,
%                   including tolerances for SVD truncation
%   - Mdom        : Global mass matrix for the mechanical domain
%   - PhiDEFbs    : Basis of basic deformation modes (e.g., affine/stretching)
%   - PhiRB       : Basis for rigid body modes (translations/rotations)
%   - MdomCHOL    : Cholesky factor of the mass matrix Mdom (such that Mdom = MdomCHOL' * MdomCHOL)
%
% OUTPUTS:
%   - UpsilonDEF     : Mass-orthonormal basis of deformation modes (including basic and nonlinear modes)
%   - S_upsilonDEF   : Diagonal matrix of singular values associated with the retained modes
%   - V_upsilonDEF   : Right singular vectors from the partitioned SVD decomposition
%
% FUNCTIONAL STEPS:
%   1. Rigid body modes are removed from each snapshot in SNAPcompl using
%      the projection operator defined by PurgeRigidBodyDisp.
%
%   2. The snapshots are mass-normalized using MdomCHOL and assembled into a
%      cell array A, where A{1} corresponds to the basic deformation modes
%      and A{2:end} to the purged nonlinear complements.
%
%   3. Each block in A is normalized with respect to the Frobenius norm to
%      ensure uniform scaling across partitions during SVD.
%
%   4. A partitioned SVD is performed on the blocks in A using the Sequential
%      Randomized SVD (SRSVD) algorithm. The SRSVD processes each block
%      sequentially, enriching the basis with left singular vectors from the
%      orthogonal complement of previous blocks. This method is particularly
%      advantageous for large-scale problems where matrix sizes approach
%      memory limitations, as it avoids the need to assemble the entire
%      snapshot matrix at once. For more details, see:
%      J.A. Hernández et al., "CECM: A continuous empirical cubature method
%      with application to the dimensional hyperreduction of parameterized
%      finite element models," arXiv:2308.03877.
%
%   5. The resulting left singular vectors are mapped back to the original
%      displacement space by inverting the mass normalization using MdomCHOL.
%
% NOTE:
%   This approach is suitable for problems involving large nonlinearities
%   and aims to extract low-rank approximations that capture both
%   fundamental and enriched strain deformation behaviors. The output
%   UpsilonDEF can be used for constructing reduced-order models in
%   multiscale or nonlinear contexts.
%   Joaquin A. Hernandez, 17-May-2025, Balmes 185, Barcelona 
%--------------------------------------------------------------------------
if nargin == 0
    load('tmp1.mat')
end

% Purging rigid-body modes
SNAPcompl = PurgeRigidBodyDisp(DATAoffline,SNAPcompl,PhiRB,Mdom,MESH) ;

DATAoffline = DefaultField(DATAoffline,'IndexesCriticalSnapshots',{});
DATAoffline = DefaultField(DATAoffline,'SVD_TOLERANCE_CriticalSnapshots',0);


if isempty(DATAoffline.IndexesCriticalSnapshots)
    
    % Partitioned SVD *********
    % -------------------------
    nproj = 1 + length(SNAPcompl) ;
    A  = cell(1,nproj) ;
    % First block = basic modes
    A{1} = MdomCHOL*PhiDEFbs ;   
    A{1} = A{1}/norm(A{1},'fro') ;  % Normalization
 
    for iproj = 1:length(SNAPcompl)
        A{iproj+1} =  MdomCHOL*SNAPcompl{iproj} ; 
         A{iproj+1}= A{iproj+1}/norm(A{iproj+1},'fro')  ;
   
    end
 
    TOL = DATAoffline.TOLSVD_complementary_modes*ones(size(SNAPcompl)) ;
    RELTOL = [0,TOL] ;
 
    DATAlocS = [] ;
    [U,S_upsilonDEF,V_upsilonDEF] = SRSVD(A,RELTOL,DATAlocS) ;
    
    UpsilonDEF = MdomCHOL\U ;  % All modes (including basic ones)
    
else
    % APPROACH INTRODUCED IN MAY 17th 2025, see
    % /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/12_REPETITIVE_TRAIN/ONE_CELL_COMPRESSION.mlx
    % Some critical snapshots are assigned TOLERANCE = 0 
    [UpsilonDEF,S_upsilonDEF,V_upsilonDEF] =  CriticalStrainModesSVD_EIFEM(SNAPcompl,DATAoffline,MdomCHOL,PhiDEFbs) ; 
    
    
    
    
end