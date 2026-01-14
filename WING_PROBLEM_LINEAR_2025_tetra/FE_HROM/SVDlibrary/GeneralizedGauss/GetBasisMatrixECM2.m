function [Lambda,VSinv] = GetBasisMatrixECM2(A,W,DATA)
% SVD of A*sqrt(W)
% --------------------------------------
%--------------------------------------------------------------------------
% function [Lambda, VSinv] = GetBasisMatrixECM2(A, W, DATA)
%
% PURPOSE:
%   Computes an orthonormal basis matrix Lambda for the column space of 
%   the weighted sampling matrix A (i.e., A * sqrt(W)) using either a standard
%   or a randomized partitioned Singular Value Decomposition (SVD), depending 
%   on the problem size and settings in DATA. The method is a crucial first 
%   step in the Empirical Cubature Method (ECM) and its continuous extension (CECM),
%   as it transforms the original integration problem into one over orthogonal basis functions.
%
% INPUTS:
%   - A     : Either a matrix [M x N] or a cell array of matrices representing 
%             the integrand functions evaluated at Gauss points.
%   - W     : [M x 1] Vector of FE quadrature weights (Gauss weights * Jacobian)
%   - DATA  : Structure with ECM settings and flags, including:
%              * TOL_SVD_A: Tolerance for SVD truncation
%              * UsePartitionedRandomizedAlgorithm: 1 to use randomized SVD
%              * Method_Evaluation_Basis_Integrand: determines behavior with constant modes
%
% OUTPUTS:
%   - Lambda : [M x r] Matrix of r orthonormal basis vectors (left singular vectors of A*sqrt(W))
%   - VSinv  : [r x N] Matrix such that A ≈ Lambda * S * V', and VSinv = inv(S) * V'
%              (used for fast projection/evaluation of integrand functions)
%
% METHOD:
%   1. Scales A with sqrt(W) to define the weighted SVD problem.
%   2. Uses one of the following:
%      - RSVDTrowblock: randomized SVD over (possibly partitioned) matrix A
%      - SVDT: classical SVD with relative truncation
%   3. Constructs the inverse projection matrix VSinv from singular values and vectors.
%   4. (Optional) Appends an additional column to Lambda to ensure exact integration 
%      of the constant function, if not already captured in the basis (important for 
%      physical conservation laws and to avoid ill-posedness).
%
% REMARKS:
%   - The function can handle both in-memory and file-based partitions of A via cell arrays.
%   - In Method_Evaluation_Basis_Integrand = 2, the function enforces exact integration 
%     of the constant function and throws an error if it is not satisfied.
%   - When Method_Evaluation_Basis_Integrand ≠ 2, the constant mode is appended only 
%     if the residual between sqrt(W) and its projection onto the computed basis 
%     exceeds a small tolerance (default TOL = 1e-10).
%
% EXAMPLE USAGE:
%   [Lambda, VSinv] = GetBasisMatrixECM2(A, wFE, struct('TOL_SVD_A', 1e-8, ...
%                                                    'UsePartitionedRandomizedAlgorithm', 1, ...
%                                                    'Method_Evaluation_Basis_Integrand', 1));
%
% SEE ALSO:
%   - ContinuousECMgen2023
%   - RSVDTrowblock (randomized SVD over row blocks)
%   - SVDT (standard truncated SVD)
%   - Hernández et al. (2024), *CECM: A continuous empirical cubature method*
%     in *Comput. Methods Appl. Mech. Engrg.*, Vol. 418, 116552
%
% AUTHOR:
%   J.A. Hernández (UPC/CIMNE), with updates for CECM framework, April–May 2024.
% Comments by ChatGPT4, 29-May-2025
%--------------------------------------------------------------------------

%dbstop('5')
if nargin == 0
    load('tmp.mat')
    
end


%
% % Exact integral
% % ---------------
% if ~iscell(A)
%     INTexac = A'*W ; % High-fidelity integral of all snapshots
%     A = bsxfun(@times,A,sqrt(W)) ;
% else
%     disp('Computing exact integrals, and multiplying by sqrt(W)')
%     INTexac = cell(size(A)) ;
%     for imat = 1:length(A)
%         disp(['iblock = ',num2str(imat)])
%         if isnumeric(A{imat})
%             INTexac{imat} =A{imat}'*W ;
%             A{imat} = bsxfun(@times,A{imat},sqrt(W)) ;
%         else
%             disp(['Retrieving from memory ...'])
%             SSS = load(A{imat}) ;
%             disp([' ... Done'])
%             fff = fieldnames(SSS) ;
%             Ai = SSS.(fff{1}) ; SSS = [] ;
%             INTexac{imat} =Ai'*W ;
%             Ai = bsxfun(@times,Ai,sqrt(W)) ;
%             disp(['Saving again in  memory ... (multiplied by sqrt(W))'])
%             save(A{imat},'Ai') ;
%             disp([' ... Done'])
%         end
%     end
%     INTexac = cell2mat(INTexac') ;
% end


% if ~iscell(A)
%     A = bsxfun(@times,A,sqrt(W)) ;
% else
%     for ipart = 1:length(A)
%         if isnumeric(A{ipart})
%             A{ipart} = bsxfun(@times,A{ipart},sqrt(W)) ;
%         else
%             SSS = load(A{imat}) ;
%             fff = fieldnames(SSS) ;
%             %    Ai = SSS.(fff{1}) ; SSS = [] ;
%             INTexac{imat} =SSS.(fff{1})'*W ;
%         end
%     end
% end


% Orthogonal basis matrix Lambda for the column space of A
disp('SVD of Aw ')
DATA = DefaultField(DATA,'UsePartitionedRandomizedAlgorithm',1) ;

% if iscell(A)
%     DATA.UsePartitionedRandomizedAlgorithm = 1 ;
%
% end
DATA = DefaultField(DATA,'TOL_SVD_A',0) ;

if  ~iscell(A)
    if DATA.UsePartitionedRandomizedAlgorithm == 1
        DATAlocSVD.HIDE_OUTPUT = 0 ;
        A = bsxfun(@times,A,sqrt(W)) ;
        tic
        [Lambda,S,V ] =      RSVDTrowblock({A},DATA.TOL_SVD_A,DATAlocSVD);
        toc
    else
        DATAsvd.RELATIVE_SVD =1 ;
        A = bsxfun(@times,A,sqrt(W)) ;
        [Lambda,S,V] =    SVDT(A,DATA.TOL_SVD_A,DATAsvd ) ;
    end
else
    %  EPSILON_GLO = DATA.TOL_SVD_A*ones(size(A)) ;
    DATAlocSVD.WeightsPremultipy_matrix = sqrt(W) ;
    
    DATAlocSVD.HIDE_OUTPUT = 0 ;
    tic
    [Lambda,S,V ] =      RSVDTrowblock(A,DATA.TOL_SVD_A,DATAlocSVD );
    toc
    disp(['Number of modes =',num2str(length(S)),' (of ',num2str(size(V,2)),'columns'])
end
VSinv = bsxfun(@times,V',1./S)' ;


disp(['Number of left singular vectors  = ',num2str(size(Lambda,2))])

% % Enlarge the basis matris for Aw
a  = sqrt(W) - Lambda*(Lambda'*sqrt(W)) ;
TOL  = 1e-10 ; TOLmax = 0.9 ;
if  DATA.Method_Evaluation_Basis_Integrand == 2
    if norm(a) > TOLmax*norm(sqrt(W))
        error('Include in the list of integrand functions the constant functions (otherwise the problem will be  ill-posed)')
    end
    
else
    if norm(a) > TOL
        INCLUDE_ADDITIONAL_COLUMN = 1;
    else
        INCLUDE_ADDITIONAL_COLUMN = 0 ;
    end
    if    INCLUDE_ADDITIONAL_COLUMN ==1
        a = a/norm(a) ;
        Lambda = [a,Lambda] ;
    end
end


