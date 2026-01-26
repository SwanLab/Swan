function [Q] = QbasisMatrixIntegrand_givenSNAPr(SNAPredFINT_nw,DATA,wSTs,DATAoffline)
%--------------------------------------------------------------------------
% QbasisMatrixIntegrand_givenSNAPr
%
% Purpose
% -------
% Build an orthonormal basis Q for the column space of *weighted* integrand
% snapshots used by ECM, with the option to treat INTERNAL FORCES and
% REACTIONS using different compression tolerances.
%
% Key idea: pre-multiply each integrand by sqrt(wSTs) and compute a reduced
% basis via (randomized) SVD so that ECM later operates on a compact
% integrand subspace. The constant mode (sqrt(wSTs)) is enforced by
% appending it if not already captured by Q (to integrate constants exactly).
%
% Signature
% ---------
%   function [Q] = QbasisMatrixIntegrand_givenSNAPr(SNAPredFINT_nw,DATA,wSTs,DATAoffline)
%
% Inputs
% ------
%   SNAPredFINT_nw : Snapshot matrix (or cell array convertible to a matrix)
%                    of *unweighted* integrands. After internal conversion
%                    this is treated as a 2-D matrix:
%                        size(SNAPredFINT_nw) = [nGP , nCOL]
%                    where nGP = number of integration points (rows after
%                    weighting by sqrt(wSTs)) and nCOL = total number of
%                    snapshot columns (possibly including both internal
%                    force components and reaction components, stacked).
%
%   DATA            : Struct with auxiliary fields:
%       .Matrix_Include_SVD   (optional, nGP×k0)  Extra columns to be
%                           forcibly included in the SVD space (e.g. known
%                           important modes). These are also weighted by
%                           sqrt(wSTs) before SVD.
%       .IndexMatrixECM.fint  (vector of column indices *within one snapshot
%                           block*) identifying internal-force columns.
%       .IndexMatrixECM.react (vector of column indices *within one snapshot
%                           block*) identifying reaction columns.
%                           When separate tolerances are requested, these
%                           two sets are expanded to all snapshots and used
%                           to split the global column set into
%                           {FINT, REACT}.
%
%   wSTs            : (nGP×1) vector of Gauss weights at the integration
%                     points. The code uses sqrt_wST = sqrt(wSTs).
%
%   DATAoffline     : Struct controlling offline compression:
%       .errorFINT                    (scalar ≥ 0) tolerance for INTERNAL
%                                     FORCES when unified tolerances are
%                                     used, or for the FINT block when
%                                     separate tolerances are used.
%       .errorFINT_reactions          (optional, scalar ≥ 0). If present,
%                                     REACTIONS are compressed with this
%                                     tolerance, enabling separate control
%                                     of accuracy for reactions.
%       .Exponent_Function_relating_global_local_TOL_fint
%                                     (optional, scalar). If non-empty,
%                                     activates the *exponential-law*
%                                     tolerance strategy (global SVD
%                                     driven by per-row exponents). This
%                                     option is currently **incompatible**
%                                     with separate reaction tolerances.
%
% Output
% ------
%   Q : (nGP×r) matrix with orthonormal columns spanning the weighted
%       snapshot space per the requested tolerances. If the constant vector
%       sqrt(wSTs) is not in span(Q), it is appended and re-normalized so
%       that constant integrands are integrated exactly by ECM.
%
% Algorithm (high level)
% ----------------------
% 1) Weighting:
%       sqrt_wST = sqrt(wSTs);
%       SNAP_w   = bsxfun(@times, SNAPredFINT_nw, sqrt_wST);
%       If Matrix_Include_SVD is provided, weight it as well.
%
% 2) SVD compression:
%    • Unified tolerance (default): run (randomized) SVD with tolerance
%      DATAoffline.errorFINT over all columns (plus any “include” columns).
%    • Split tolerances (internal forces + reactions): build global column
%      index sets for FINT and REACT by repeating the per-block indices
%      over the snapshot blocks; run a multi-block SVD with tolerances
%      [errorFINT, errorFINT_reactions].
%    • Exponential-law option: if
%      Exponent_Function_relating_global_local_TOL_fint is non-empty,
%      call Qmatrix_for_ECM_exponentialLAWtolerances. (Not allowed together
%      with .errorFINT_reactions.)
%
% 3) Constant-mode enforcement:
%       a = sqrt_wST - Q*(Q'*sqrt_wST);
%       if norm(a) > 1e-10, set Q = [Q, a/norm(a)].
%
% Notes & Conventions
% -------------------
% • Dimensions: rows (nGP) index quadrature points; columns index snapshot
%   realizations and/or components. The function assumes column grouping by
%   “snapshot blocks,” so that DATA.IndexMatrixECM.fint/react define the
%   within-block columns and are replicated across snapshots internally.
%
% • Tolerances: setting .errorFINT = 0 disables truncation for the FINT
%   block (a diagnostic SVD-error plot is available but off by default).
%
% • Include-matrix: columns in DATA.Matrix_Include_SVD are *prepended* to
%   the SVD with zero tolerance in the multi-block call, i.e., they are
%   always kept before applying the tolerances to the snapshot blocks.
%
% • Separate tolerances are useful when reaction magnitudes/scales differ
%   from internal forces, preventing one block from dominating the SVD and
%   allowing tighter control of each error budget.
%
% • Error handling: requesting both the exponential-law option and
%   .errorFINT_reactions triggers an explicit error, since the current
%   implementation only supports split tolerances under the standard SVD
%   path.
%
% Typical usage
% -------------
%   DATA.IndexMatrixECM.fint  = [1 3 5];      % columns per block (example)
%   DATA.IndexMatrixECM.react = [2 4];        % columns per block (example)
%   DATAoffline.errorFINT           = 1e-6;
%   DATAoffline.errorFINT_reactions = 1e-4;   % enable split tolerances
%   Q = QbasisMatrixIntegrand_givenSNAPr(SNAP, DATA, wSTs, DATAoffline);
%
% Author / History
% ----------------
%   Adapted from QbasisMatrixIntegrand_givenSNAP to support split
%   tolerances for internal forces vs. reactions.
%   J.A. Hernández (UPC), 9-Oct-2025, Barcelona (Balmes 185).
%--------------------------------------------------------------------------


if nargin == 0
    load('tmp1.mat')
    %     DATA.Matrix_Include_SVD = [] ;
    %     %  DATAoffline.Exponent_Function_relating_global_local_TOL_fint =[];
    %     % DATAoffline.errorFINT = 1e-6 ;
    %     DATAoffline.errorFINT_reactions = 1e-4;
end

%SNAPredFINT_nw = BasisF_from_BasisStress_PK1(BstRED_l,BasisPone,DATA)  ;
%  wSTs_LOC = OPERFE.wSTs ;
sqrt_wST = sqrt(wSTs) ;

DATAoffline = DefaultField(DATAoffline,'Exponent_Function_relating_global_local_TOL_fint',0) ;% = 0.01;
DATA = DefaultField(DATA,'Matrix_Include_SVD',[]) ;% = 0.01;


DATAoffline = DefaultField(DATAoffline,'errorFINT_reactions',[]) ;


info = whos('SNAPredFINT_nw');
size_in_bytes = info.bytes;
size_in_MB = size_in_bytes / (1024^2);
fprintf('SNAPredFINT_nw uses %.2f MB of memory.\n', size_in_MB);
DATAoffline = DefaultField(DATAoffline,'LimitSizeMBytesA_fint',1000) ;


if size_in_MB > DATAoffline.LimitSizeMBytesA_fint
    % Cell array  SNAPredFINT_nw should be treated as cell array
    Q = Qmatrix_LOC_forECM_cell(SNAPredFINT_nw,DATA,wSTs,DATAoffline,sqrt_wST) ;
    
else
    
    Q = Qmatrix_LOC_forECM(SNAPredFINT_nw,DATA,wSTs,DATAoffline,sqrt_wST) ;
end

%end
if DATAoffline.errorFINT == 0
    ifig = 3000 ;
    plot_svd_error = 0 ;
    if plot_svd_error == 1
        SVDplotERROR_local(S,ifig) ;
    end
end

% % Enlarge the basis matris for SNAPredFINT
a  = sqrt_wST - Q*(Q'*sqrt_wST) ;
if norm(a) > 1e-10
    a = a/norm(a) ;
    Q = [Q,a] ;
    
    %     aNW = a./sqrt_wST ;
    %     SNAPredFINT_nw = [aNW,SNAPredFINT_nw] ; %
    
end