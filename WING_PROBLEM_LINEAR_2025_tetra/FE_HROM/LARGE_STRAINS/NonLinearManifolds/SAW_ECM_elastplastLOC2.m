function [ECMdata] = ...
    SAW_ECM_elastplastLOC2(SNAPfint_lin,SNAPfint_non,DATA,wSTs,DATAoffline,qPLAST)
%--------------------------------------------------------------------------
%==========================================================================
% SAW_ECM_elastplastLOC2
%--------------------------------------------------------------------------
% Purpose
%   Subspace–Adaptive Weights Empirical Cubature (SAW-ECM) point candidate
%   generator tailored for elastoplastic (nonlinear) internal forces.
%   Given elastic and nonlinear snapshots of internal forces, it:
%     1) builds a weighted elastic basis (plus a constant mode),
%     2) for each plastic “cluster” (time step / load step), extracts a
%        local weighted basis that is orthogonal to the elastic subspace
%        while allowing neighborhood overlap,
%     3) passes those per-cluster local bases to SAW_ECM to obtain the
%        candidate quadrature points used later by an ECM/DECM solver.
%
% Inputs
%   SNAPfint_lin : cell array {nBlocks} of matrices (nGP × nSnap_e,i)
%                  Elastic (linearized) internal-force snapshots per block
%                  (or concatenate-able batches). Each column is a snapshot.
%   SNAPfint_non : cell array {nClusters} of matrices (nGP × nSnap_n,i)
%                  Nonlinear (elasto-plastic) internal-force snapshots for
%                  each cluster (e.g., time/load step neighborhood).
%   DATA         : struct with mesh/problem info (passed through to SAW_ECM).
%   wSTs         : vector (nGP × 1) of positive weights for the stress
%                  integration points (Gauss weights × detJ × any scalings).
%   DATAoffline  : struct with offline settings/flags. Recognized fields:
%       .ECM_Number_Snapshots_Overlapping  (int, default = 1)
%            Number of neighboring clusters to include on each side when
%            extracting the local basis of a cluster (overlapping window).
%       .IncludeElastic_SAW_ECM            (bool, default = 1)
%            If true, uses elastic snapshots to build elastic subspace.
%       .IncludeConstantFunction_SAW_ECM   (bool, default = 1)
%            If true, appends a constant function to the elastic subspace.
%       .errorFINT                         (double, default = 1e-4)
%            Tolerance for capturing cluster-internal (nonlinear) modes.
%       .errorFINT_neightCLUSTER           (double, default = 1e-3)
%            Tolerance for capturing neighboring-cluster contributions.
%       (Other fields can be present; they are forwarded to SAW_ECM.)
%   qPLAST       : (optional) plastic variables/state info (unused here, kept
%                  for interface compatibility / future extensions).
%
% Outputs
%   ECMdata_cluster   : struct returned by SAW_ECM containing per-cluster
%                       metadata for empirical cubature (indices, weights,
%                       errors, timings, etc. depending on SAW_ECM impl.).
%   setCandidates     : vector of unique candidate Gauss-point indices (1-based)
%                       selected across all clusters (union).
%   LIST_OF_CANDIDATES: cell array {nClusters}; each entry lists the candidate
%                       indices proposed for that cluster before taking unions.
%
% Key steps (high level)
%   • Weighting:
%       Wfe = diag(wSTs). The function uses weighted SVD operators WSVDT/SRSVD.
%       Internally, WSVDT returns a Cholesky-like factor sqW such that
%       sqW'*sqW ≈ Wfe (used to keep weighting consistent across calls).
%
%   • Elastic subspace:
%       - Concatenate elastic snapshots: SNAPfint_lin = cell2mat(SNAPfint_lin).
%       - Compute truncated SVD via SRSVD(SNAPfint_lin, TOL) with TOL = 1e-6.
%       - Apply weighted orthonormalization WSVDT(·, Wfe) to get PhiELAS.
%       - Optionally append a constant function (all ones) and re-orthonormalize
%         with WSVDT([PhiELAS, ones(...)], [], DATAloccc) using the stored sqW.
%
%   • Nonlinear local (cluster) subspaces with overlap:
%       For each cluster i:
%         - Build the overlapping window [i-k, …, i, …, i+k], where
%           k = DATAoffline.ECM_Number_Snapshots_Overlapping, clipped to bounds.
%         - Split window into “current cluster” SNAPlocI and its neighbors
%           SNAPlocNEIGTH.
%         - Compute a weighted SRSVD on the stacked set
%             { sqW*PhiELAS, sqW*SNAPlocI, sqW*SNAPlocNEIGTH }
%           with per-block tolerances
%             [ 0,  TOL_icluster,  TOL_neigh_All ]
%           where TOL_icluster = DATAoffline.errorFINT and
%                 TOL_neigh_All = DATAoffline.errorFINT_neightCLUSTER (vector).
%         - Store the resulting left singular vectors UU (weighted orthonormal
%           basis) as the local basis SNAPfint_w{i}.
%         - Track NumberModesCluster(i) = length(SS) for diagnostics/plots.
%
%   • Candidate selection:
%       - Call SAW_ECM(wSTs, DATAoffline, SNAPfint_w, DATA) to compute per-cluster
%         candidate sets and the global union (setCandidates).
%
% External dependencies (expected on path)
%   DefaultField  : utility to set default struct fields if missing.
%   SRSVD         : (block-)tolerant SVD with optional per-set truncation.
%   WSVDT         : weighted SVD / weighted orthonormalization (returns sqW).
%   SAW_ECM       : subspace-adaptive weights empirical cubature driver that
%                   produces candidate Gauss points from per-cluster bases.
%
% Notes
%   • nargin==0 branch loads a demo dataset (tmp2.mat) and sets convenient
%     default tolerances and flags for quick experimentation.
%   • The plotting (figure 42) is purely diagnostic to visualize the number
%     of retained internal-force modes per cluster.
%   • The function focuses on *candidate generation*. Actual reduced
%     quadrature weights and final selection for integration are computed
%     inside SAW_ECM / downstream ECM routines.
%   • Overlap level controls robustness vs. locality: higher overlap better
%     preserves trajectories that leak across clusters, at increased cost.
%   • Tolerances:
%       - errorFINT (cluster) should be tighter than errorFINT_neightCLUSTER
%         (neighbors) to prioritize fidelity on the current cluster while
%         still capturing essential neighbor content.
%   • Constant mode:
%       Including a constant internal-force mode helps preserve equilibrium
%       offsets and improves stability in subsequent ECM selection.
%
% Complexity (rough guide)
%   Let nGP be the number of Gauss points, and r be typical retained rank.
%   Each SRSVD over the stacked sets costs ~O(nGP·r^2) per cluster (assuming
%   truncated algorithms). Overlap increases stacked size but improves locality.
%
% Diagnostics printed to console
%   - Number of elastic modes retained
%   - Number of candidate points and their indices (global union)
%
% Example (pseudo-usage)
%   DATAoffline.ECM_Number_Snapshots_Overlapping = 1;
%   DATAoffline.IncludeElastic_SAW_ECM          = 1;
%   DATAoffline.IncludeConstantFunction_SAW_ECM = 1;
%   DATAoffline.errorFINT                       = 1e-4;
%   DATAoffline.errorFINT_neightCLUSTER         = 1e-3;
%   [ECMdata_cluster,setCandidates,LIST] = SAW_ECM_elastplastLOC2( ...
%       SNAPfint_lin, SNAPfint_non, DATA, wSTs, DATAoffline, qPLAST);
%
% References (context)
%   • Empirical Cubature / Hyper-reduction literature by Hernández et al.,
%     Bravo et al., de Parga, Rossi, et al. (SAW-ECM, DECM/ECM variants).
%   • Weighted SVD / orthonormalization for integration-consistent bases.
%
% Version
%   This variant includes neighborhood overlap for nonlinear clusters and
%   appends a constant function to the elastic basis when requested.
%==========================================================================

%--------------------------------------------------------------------------
if nargin == 0
    load('tmp3.mat')
    %     DATAoffline.ECM_Ratio_Number_Clusters_Snapshots = 1 ;
    DATAoffline.ECM_Number_Snapshots_Overlapping =20 ;
    
    DATAoffline.errorFINT  =1e-3;
 %   DATAoffline.errorFINT_neightCLUSTER = 1e-3;
   close all
end


DATAoffline = DefaultField(DATAoffline,'ECM_Number_Snapshots_Overlapping',1) ;
  
SNAPfint_lin = cell2mat(SNAPfint_lin) ;
TOL = 1e-6;
[PhiELAS,S,V] = SRSVD(SNAPfint_lin,TOL) ;
%Wfe = diag(sparse(wSTs)) ;
Wfe = diag(sparse(wSTs)) ;
[PhiELAS,ssss,~,sqW] = WSVDT(PhiELAS,Wfe) ;
disp(['Number of elastic modes (internal forces) = ',num2str(size(PhiELAS,2))]);
%  Enlarged to include the constant function
DATAloccc.Mchol = sqW ;
%
% DATAoffline = DefaultField(DATAoffline,'IncludeConstantFunction_SAW_ECM',1)  ;
%
% if DATAoffline.IncludeConstantFunction_SAW_ECM == 1
PhiELAS = WSVDT([PhiELAS,ones(size(PhiELAS,1),1) ],[],DATAloccc ) ;
%end
%%%% Next we compute the orthogonal complement
SNAPfint_w = SNAPfint_non;
TOL_icluster = DATAoffline.errorFINT ;
NumberModesCluster = zeros(size(SNAPfint_non)) ;
DATADD.HIDE_OUTPUT = 1;

%DATAoffline.errorFINT  =1e-4;
DATAoffline= DefaultField(DATAoffline,'errorFINT_neightCLUSTER',DATAoffline.errorFINT) ;
TOL_neigh  = DATAoffline.errorFINT_neightCLUSTER ;
DATADD.SortByTolerances = 0 ;

% ISSUE WITH TOLERANCES 
% WE WISH APPROXIMATIONS TO REFLECT THE IMPORTANCE OF PLASTIC SNAPSHOTS
% ELASTIC MODES MUST BE EXACTLY CAPTURED. THEY ARE NORMALIZED SUCH THAT
% PhiELAST'*W*PhiELAST = 1 
% What'we are going to do is to compute the same norm for the last plastic
% snapshots
normLAST = norm(SNAPfint_non{end}'*Wfe*SNAPfint_non{end},'fro')/sqrt(size(SNAPfint_non{end},2)) ;
 


for icluster = 1:length(SNAPfint_non)
 
    disp(['icluster = ',num2str(icluster)])
    icluster_back = max(1,icluster-DATAoffline.ECM_Number_Snapshots_Overlapping) ;
    icluster_forw = min(length(SNAPfint_non),icluster+DATAoffline.ECM_Number_Snapshots_Overlapping) ;
    icluster_select = [icluster_back:icluster_forw] ;
    
    [indNEIGH_loc,indNEIGH_glo] = setdiff(icluster_select,icluster) ;
   % TOL_neigh_All = ones(size(indNEIGH_glo))*TOL_neigh ;
    % TOL_plast(indNEIGH_loc) = TOL_neigh ;
    %   [indCLUST_loc,indCLUST_glo] = setdiff(icluster_select,indNEIGH_glo) ;
    
        SNAPlocI = cell2mat(SNAPfint_non(icluster)) ;

    
   % normLOC = norm(SNAPlocI'*Wfe*SNAPlocI,'fro')/sqrt(size(SNAPlocI,2)) ;
%     factor_tolerance = normLAST/normLOC ; 
%     TOL_cluster_weigthed = min(factor_tolerance*TOL_icluster,0.99) ; 
%     TOL_neigh_we = min(factor_tolerance*TOL_neigh,0.99) ; 
     
    % DATADD.EPSILON_ABSOLUTE = 1; 
    %[UU,SS,~] =   SRSVD({sqW*PhiELAS,sqW*SNAPlocI,sqW*SNAPlocNEIGTH},[0,TOL_cluster_weigthed,TOL_neigh_we],DATADD);
    
%     SNAPlocI = SNAPlocI/norm(SNAPlocI,'fro') ; 
     SNAPlocNEIGTH = cell2mat(SNAPfint_non(indNEIGH_glo))  ; 
%     SNAPlocNEIGTH = SNAPlocNEIGTH/norm(SNAPlocNEIGTH,'fro') ; 
    
    PhiPlast = [SNAPlocI,SNAPlocNEIGTH] ; 
    PhiPlast = PhiPlast -PhiELAS*(PhiELAS'*Wfe*PhiPlast) ; 
    
    DATADDDD.ISRELATIVE = 1; 
   % PhiELAS = 10*PhiELAS/norm(PhiELAS,'fro');
    [UU,SS,VV] =  SVDT(PhiPlast,DATAoffline.errorFINT,DATADDDD); 
    UUf =  SVDT([sqW*PhiELAS,sqW*UU]); 
    
    
    SNAPfint_w{icluster} =  UUf ;
    NumberModesCluster(icluster) = length(SS) ;
end

figure(42)
hold on
xlabel('Snapshot')
ylabel('Number of internal force modes')
bar(NumberModesCluster)
%
% else
%      TOL = DATAoffline.errorFINT ;
%        NumberModesCluster = zeros(size(SNAPfint_non)) ;
%       for isnap = 1:length(SNAPfint_non)
%         [UU,SS,~] =   SRSVD(sqW*SNAPfint_non{isnap},TOL);
%         SNAPfint_w{isnap} =  UU ;
%         NumberModesCluster(isnap) = length(SS) ;
%     end
%
%     figure(42)
%     hold on
%     xlabel('Snapshot')
%     ylabel('Number of internal force modes')
%     bar(NumberModesCluster)
% end

%
[ECMdata_cluster,setCandidates,LIST_OF_CANDIDATES] = ...
    SAW_ECM(wSTs,DATAoffline,SNAPfint_w,DATA) ;

disp(['Number of Candidate Points'])
disp(length(setCandidates))
disp(['setCandidates = '])
disp(setCandidates')


%%%%%%%%%%%%%%%%%5
%Creating the chart qPLAST-wECM,zECM
setPointsALL = cell(size(ECMdata_cluster)) ;

for icluster = 1:length(ECMdata_cluster)
    setPointsALL{icluster} =  ECMdata_cluster{icluster}.setPoints  ;
end
setPointsALL  =unique(cell2mat(setPointsALL' )) ;

ncluster = length(ECMdata_cluster) ;
npointsALL = length(setPointsALL) ;
wALL = zeros(npointsALL,ncluster) ;

for icluster = 1:ncluster
    setPloc = ECMdata_cluster{icluster}.setPoints  ;
    [dummy1,III,JJJ] = intersect(setPloc,setPointsALL,'stable') ;
    wALL(JJJ,icluster) =  ECMdata_cluster{icluster}.wRED  ;
end

ECMdata.setPoints = setPointsALL ;
ECMdata.wRED.Values = wALL ;
ECMdata.wRED.q = qPLAST ;
setElements_cand = large2smallINCLUDEREP(setPointsALL,DATA.MESH.ngaus) ;
ECMdata.setElements = setElements_cand ;

ECMdata.wRED.IndexDOFl_q = 2;

% 
% warning('BORRAR ESTO, ESTOY ENGAÑANDO AL PROGRAMA')
%     load('tmpWORKecm.mat','setPoints','wRED')
% ECMdata.setPoints = setPoints ;
% ECMdata.wRED.Values = repmat(wRED,1,length(qPLAST)) ;