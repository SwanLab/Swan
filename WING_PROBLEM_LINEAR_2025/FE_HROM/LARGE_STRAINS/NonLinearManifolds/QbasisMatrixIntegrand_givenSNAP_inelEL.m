function [Q,SNAPredFINT_nw] = QbasisMatrixIntegrand_givenSNAP_inelEL(SNAPredFINT_nw_lin,SNAPredFINT_nw_non ,DATA,wSTs,DATAoffline)
%--------------------------------------------------------------------------
% This is a modification of QbasisMatrixIntegrand_givenSNAP (described below). Now we have
% two distinct contributions:  one linear, and one nonlinear
% The linear one is to be exactly integrated
% JAHO, 25th August 2025, Balmes 185, Barcelona.
%
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
    load('tmp.mat')
    % DATAoffline.errorFINT = 1e-6 ;
end

%SNAPredFINT_nw = BasisF_from_BasisStress_PK1(BstRED_l,BasisPone,DATA)  ;
%  wSTs_LOC = OPERFE.wSTs ;
sqrt_wST = sqrt(wSTs) ;


DATAoffline = DefaultField(DATAoffline,'IncludeLinearSnapshots_ECM',1) ;

if  DATAoffline.IncludeLinearSnapshots_ECM == 1
   
%     SNAPredFINT_nw = {SNAPredFINT_nw_lin,SNAPredFINT_nw_non} ;
%     TOL_loc = [0,DATAoffline.errorFINT] ;
    
     SNAPredFINT_nw = {SNAPredFINT_nw_non,SNAPredFINT_nw_lin} ;
    TOL_loc = [DATAoffline.errorFINT,0] ;
    
else
     warning('This option has proved to be not robust')
    SNAPredFINT_nw = {SNAPredFINT_nw_non} ;
    TOL_loc = [DATAoffline.errorFINT] ;
end

SNAPredFINT = cell(size(SNAPredFINT_nw));
for iproj = 1:length(SNAPredFINT_nw)
    SNAPredFINT{iproj} = bsxfun(@times,SNAPredFINT_nw{iproj},sqrt_wST) ;
end

DATAsvdLOC.HIDE_OUTPUT =  1;
DATAsvdLOC.SortByTolerances = 0; 
[Q,S,V] = SRSVD(SNAPredFINT,TOL_loc,DATAsvdLOC) ;



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