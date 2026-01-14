function [Q] = QbasisMatrixIntegrand_givenSNAP(SNAPredFINT_nw,DATA,wSTs,DATAoffline)
%--------------------------------------------------------------------------
% function [Q, SNAPredFINT_nw] = QbasisMatrixIntegrand_givenSNAP(SNAPredFINT_nw, DATA, wSTs, DATAoffline)
%
% Purpose:
%   Constructs an orthonormal basis matrix Q for the column space of the
%   weighted internal force snapshots, used in the Empirical Cubature Method (ECM).
%
%   The function applies a square-root weighting to the integrand (internal force snapshots)
%   and computes a reduced orthonormal basis using a randomized SVD (SRSVD).
%   Optionally, the constant vector (vector of square root weights) is appended to Q
%   to ensure exact integration of constant terms.
%
% Inputs:
%   - SNAPredFINT_nw : Cell array of internal force snapshot matrices (unweighted)
%   - DATA           : Simulation and mesh configuration
%   - wSTs           : Vector of Gauss integration weights
%   - DATAoffline    : Offline configuration structure, containing:
%       • .errorFINT : tolerance used to compress the integrand via SVD
%
% Outputs:
%   - Q              : Orthonormal basis matrix for the weighted internal force snapshots
%   - SNAPredFINT_nw : Possibly augmented snapshot matrix (if constant mode was appended)
%
% Procedure:
%   1. Apply square-root weighting to each snapshot matrix.
%   2. Perform randomized SVD (SRSVD) on the weighted snapshots to obtain Q.
%   3. If needed, append a normalized residual to ensure exact integration
%      of the square-root weight vector (e.g., when Q fails to span it).
%
% Notes:
%   - If DATAoffline.errorFINT = 0, no truncation is applied.
%   - The final basis Q is guaranteed to integrate constant functions exactly,
%     provided the square-root weights lie in the column space of Q.
%   - This function is typically used in the offline stage before applying ECM.
%
% Author:
%   Joaquín A. Hernández, UPC, 12-JUL-2025
%
%--------------------------------------------------------------------------

if nargin == 0
    load('tmp1.mat')
    %DATA.Matrix_Include_SVD = [] ; 
  %  DATAoffline.Exponent_Function_relating_global_local_TOL_fint =[];
    % DATAoffline.errorFINT = 1e-6 ;
end

%SNAPredFINT_nw = BasisF_from_BasisStress_PK1(BstRED_l,BasisPone,DATA)  ;
%  wSTs_LOC = OPERFE.wSTs ;
sqrt_wST = sqrt(wSTs) ;

DATAoffline = DefaultField(DATAoffline,'Exponent_Function_relating_global_local_TOL_fint',0) ;% = 0.01;
DATA = DefaultField(DATA,'Matrix_Include_SVD',[]) ;% = 0.01;


if isempty(DATAoffline.Exponent_Function_relating_global_local_TOL_fint)
    % Global SVD
    SNAPredFINT_nw  = bsxfun(@times,cell2mat(SNAPredFINT_nw),sqrt_wST) ;
    
    
    
    
    if ~isempty(DATA.Matrix_Include_SVD)
        Matrix_Include_SVD_w  = bsxfun(@times,DATA.Matrix_Include_SVD,sqrt_wST) ;
        DATAsvdLOC.HIDE_OUTPUT =  0;
        [Q,S,V] = SRSVD({Matrix_Include_SVD_w,SNAPredFINT_nw},[0,DATAoffline.errorFINT],DATAsvdLOC) ;
    else
        
        DATAsvdLOC.HIDE_OUTPUT =  0;
        [Q,S,V] = SRSVD(SNAPredFINT_nw,DATAoffline.errorFINT,DATAsvdLOC) ;
    end
    
    
else
    
     if ~isempty(DATAoffline.errorFINT_reactions)
        error('Option only available for unified treatment of internal forces and reactions')
    end
    [Q,S,V] =  Qmatrix_for_ECM_exponentialLAWtolerances(SNAPredFINT_nw,sqrt_wST,DATAoffline) ;
end




%else
%    [Q,S,V,eSVD] = SVDT(SNAPredFINT,DATAoffline.errorFINT,DATAsvd) ;

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