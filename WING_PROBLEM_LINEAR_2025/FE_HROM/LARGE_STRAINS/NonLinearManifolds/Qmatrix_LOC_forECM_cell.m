function Q = Qmatrix_LOC_forECM_cell(SNAPredFINT_nw,DATA,wSTs,DATAoffline,sqrt_wST)
if nargin == 0
    load('tmp1.mat')
end
%--------------------------------------------------------------------------
% Qmatrix_LOC_forECM_cell
%
% Purpose
% -------
% Construct an orthonormal basis matrix Q for a collection of integrand
% snapshots provided as a **cell array**, applying weighting by the square
% root of the integration weights and optional inclusion of additional
% matrices before SVD compression.
%
% This function is a cell-oriented, memory-safe variant of
% QbasisMatrixIntegrand_givenSNAPr, designed for situations where the
% integrand snapshots are too large to concatenate into a single matrix.
% Each cell of SNAPredFINT_nw is treated as an independent block during
% randomized SVD (SRSVD) computation.
%
% Signature
% ---------
%   function Q = Qmatrix_LOC_forECM_cell(SNAPredFINT_nw, DATA, wSTs, DATAoffline, sqrt_wST)
%
% Inputs
% ------
%   SNAPredFINT_nw : cell array {nBlocks×1}
%       Each cell contains a matrix of *unweighted* integrand snapshots for
%       a subset of Gauss points or components.  The function processes
%       these blocks independently via SRSVD without concatenating them.
%
%   DATA : struct with configuration fields
%       .Matrix_Include_SVD   (optional, nGP×k0)
%           Matrix whose columns are forcibly included in the final SVD
%           subspace.  Typically used to preserve specific physical modes
%           or constraints.  If non-empty, it is automatically prepended to
%           the list of SVD blocks with zero truncation tolerance.
%       .IndexMatrixECM.fint  (optional) Indices of internal-force entries.
%       .IndexMatrixECM.react (optional) Indices of reaction entries.
%           The split-tolerance mode for (fint, react) is **not yet
%           implemented** in this version.
%
%   wSTs : vector (nGP×1)
%       Gauss integration weights at each point.  Used to compute
%       sqrt_wST if not passed explicitly.
%
%   DATAoffline : struct controlling offline compression
%       .errorFINT              : SVD truncation tolerance for the
%                                 integrand snapshots (mandatory).
%       .errorFINT_reactions    : optional tolerance for reaction blocks
%                                 (currently not implemented here).
%       .Exponent_Function_relating_global_local_TOL_fint :
%                                 optional exponential-law control
%                                 parameter. If non-empty, this version
%                                 raises an error since the exponential-law
%                                 path is not supported for cell inputs.
%
%   sqrt_wST : vector (nGP×1)
%       Precomputed square roots of the integration weights. If omitted,
%       computed internally as sqrt(wSTs).  Passed explicitly here for
%       efficiency when repeatedly called within the same offline loop.
%
% Output
% ------
%   Q : matrix (nGP×r)
%       Orthonormal basis spanning the weighted integrand space. The number
%       of columns r is determined automatically from the SVD truncation
%       tolerance.
%
% Algorithmic outline
% -------------------
% 1) Configure DATAsvdLOC structure:
%       DATAsvdLOC.HIDE_OUTPUT = 0;
%       DATAsvdLOC.WeightsPremultipy_matrix = sqrt_wST;
%
% 2) If DATA.Matrix_Include_SVD exists:
%       - prepend it to the list of snapshot blocks
%       - set its truncation tolerance to 0
%       - perform a blockwise randomized SVD using:
%             [Q,S,V] = SRSVD([{DATA.Matrix_Include_SVD}, SNAPredFINT_nw(:)'], ...
%                              DATAoffline.errorFINT, DATAsvdLOC);
%       - optionally augment Q with any residual components of
%         Matrix_Include_SVD not captured by SRSVD.
%
% 3) Otherwise:
%       - call directly
%             [Q,S,V] = SRSVD(SNAPredFINT_nw, DATAoffline.errorFINT, DATAsvdLOC);
%
% 4) If split tolerances (errorFINT / errorFINT_reactions) are provided,
%       an error is raised since this feature is not yet implemented for
%       the cell-based variant.
%
% 5) Return Q as the local orthonormal basis for use in ECM preprocessing.
%
% Notes
% -----
% • Designed for large-scale problems where concatenating all snapshot
%   blocks into a single matrix would exceed available memory.
%
% • The optional inclusion matrix is projected and re-orthogonalized if not
%   fully captured by Q, ensuring that pre-specified modes remain exactly
%   representable.
%
% • The function focuses on the unified tolerance case
%   (.errorFINT_reactions = []) and raises explicit errors for unsupported
%   options such as exponential-law or separate reaction tolerances.
%
% • The WeightsPremultipy_matrix field in DATAsvdLOC ensures that SRSVD
%   internally premultiplies each block by diag(sqrt_wST), equivalent to
%   weighting the integrand before basis extraction.
%
% Author / History
% ----------------
%   J.A. Hernández (UPC-CIMNE)
%   Barcelona, 12-Oct-2025
%   Adapted from QbasisMatrixIntegrand_givenSNAPr to handle cell-based,
%   blockwise snapshot processing with weighting and include-matrix support.
%--------------------------------------------------------------------------



%[Q] = QbasisMatrixIntegrand_givenSNAPr(SNAPredFINT_nw,DATA,wSTs,DATAoffline)
 
if isempty(DATAoffline.Exponent_Function_relating_global_local_TOL_fint)
    % Global SVD
    % SNAPredFINT_nw  = bsxfun(@times,cell2mat(SNAPredFINT_nw),sqrt_wST) ;
    if isempty(DATAoffline.errorFINT_reactions)
        
        
   %     DATA.Matrix_Include_SVD = rand(size(SNAPredFINT_nw{1})) ; warning('borrar esto')
        DATAsvdLOC.HIDE_OUTPUT =  0;
        DATAsvdLOC.WeightsPremultipy_matrix = sqrt_wST ;
        
        
        
        if ~isempty(DATA.Matrix_Include_SVD)
     
           [Q,S,V] = SRSVD( [{DATA.Matrix_Include_SVD}, SNAPredFINT_nw(:)'], DATAoffline.errorFINT, DATAsvdLOC );
           
           % % Enlarge the basis matris for SNAPredFINT
           a  = DATA.Matrix_Include_SVD - Q*(Q'*DATA.Matrix_Include_SVD) ;
           if norm(a)/norm(DATA.Matrix_Include_SVD) > 1e-10
               a = SVDT(a) ;
               Q = [Q,a] ;
                
           end
           
           
        else
            if  norm(SNAPredFINT_nw{1},'fro') <1e-10 
                SNAPredFINT_nw(1) = [] ; 
            end
            
            [Q,S,V] = SRSVD(SNAPredFINT_nw,DATAoffline.errorFINT,DATAsvdLOC) ;
        end
        
    else
        
        error('Option not implemented for this scenario....')
        
        IndLOC = DATA.IndexMatrixECM.fint ; %  IndLOC =1
        
        IndLOC = repmat(IndLOC,nsnap,1) ;
        nentriesLOC = length(DATA.IndexMatrixECM.fint) + length(DATA.IndexMatrixECM.react) ;
        FactorsM = (1:nentriesLOC:nentriesLOC*nsnap)-1 ;
        IndLOC = bsxfun(@plus,IndLOC,FactorsM');
        IndFINT = IndLOC(:) ;
        
        
        IndLOC = DATA.IndexMatrixECM.react ;  % Ind_react = 2:5
        IndLOC = bsxfun(@plus,IndLOC,FactorsM');
        IndREACT = IndLOC(:) ;
        
        
        
        
        if ~isempty(DATA.Matrix_Include_SVD)
            Matrix_Include_SVD_w  = bsxfun(@times,DATA.Matrix_Include_SVD,sqrt_wST) ;
            DATAsvdLOC.HIDE_OUTPUT =  0;
            [Q,S,V] = SRSVD({Matrix_Include_SVD_w,SNAPredFINT_nw(:,IndFINT),SNAPredFINT_nw(:,IndREACT)},...
                [0,DATAoffline.errorFINT,DATAoffline.errorFINT_reactions],DATAsvdLOC) ;
            
            % 
            
            
            
        else
            
            DATAsvdLOC.HIDE_OUTPUT =  0;
            [Q,S,V] = SRSVD({SNAPredFINT_nw(:,IndFINT),SNAPredFINT_nw(:,IndREACT)},[DATAoffline.errorFINT,DATAoffline.errorFINT_reactions],DATAsvdLOC) ;
            
            % [Q,S,V] = SRSVD([SNAPredFINT_nw,SNAPreact_nw],[DATAoffline.errorFINT],DATAsvdLOC) ;
        end
        
    end
    
else
    error('Option not implemented for such a scenario ....')
    %     if ~isempty(DATAoffline.errorFINT_reactions)
    %         error('Option only available for unified treatment of internal forces and reactions')
    %     end
    %     [Q,S,V] =  Qmatrix_for_ECM_exponentialLAWtolerances(SNAPredFINT_nw,sqrt_wST,DATAoffline) ;
    
end




%else
%    [Q,S,V,eSVD] = SVDT(SNAPredFINT,DATAoffline.errorFINT,DATAsvd) ;