function [BmatIstE,WecmE,KcLINlocE] = Bmat_weights_EIFEbubCOROT(DATA,TRANSF_COORD,EIFEoper_all)
%--------------------------------------------------------------------------
%  Bmat_weights_EIFEbubCOROT
%
%  Computes the **interscale B-matrix**, **integration weights**, and the
%  **local linear stiffness matrix** of a coarse-scale element in the 
%  co-rotational EIFEM formulation (small strains / large rotations).
%
%  This function implements the transformation of offline precomputed operators 
%  (in the reference parent domain) into the physical configuration of the
%  element, accounting for **rotation**, **uniform scaling**, and **translation**.
%
%  INPUTS:
%    - DATA          : user-defined structure with flags and configuration
%    - TRANSF_COORD  : mapping info from the reference to the physical element
%                      > includes ROTATION, SCALEFACTOR (λ), detJe, etc.
%    - EIFEoper_all  : array of precomputed parent domain operators
%
%  OUTPUTS:
%    - BmatIstE : Elemental interscale B-matrix for the current element
%    - WecmE    : Weights for the coarse-scale internal force cubature
%    - KcLINlocE: Local coarse-scale linear stiffness matrix
%
%  SCALING TRANSFORMATIONS:
%    - Interscale B-matrix:
%         BmatIstE = (1/λ) * BmatIrefST
%         (Eq. 560 in *EIFEM_largeROTfinal.pdf*):contentReference[oaicite:0]{index=0}
%    - Cubature weights:
%         WecmE = detJe * Wref
%         where detJe = λ^ndim (Eq. 561):contentReference[oaicite:1]{index=1}
%    - Linear stiffness matrix:
%         KcLINlocE = λ^(ndim-2) * KcLINref
%         (Eq. 470, related to rotated stiffness computation):contentReference[oaicite:2]{index=2}
%
%  NOTES:
%    - The parent domain index is chosen during the mapping phase
%      (`TRANSF_COORD.IndexParentDomain`).
%    - This function assumes that uniform scaling is valid
%      (`DATA.UNIFORM_SCALING_REFERENCE_ELEMENT == 1`).
%    - Rotation matrices and centroids are provided externally via `TRANSF_COORD`.
%
%  ROLE IN EIFEM:
%    These transformed matrices are used during the online assembly of:
%      - Internal force vector
%      - Linear stiffness matrix
%      - Residual compatibility condition
%
%  Author:
%    Joaquín A. Hernández Ortega, UPC/CIMNE
%    Balmes 185, Barcelona
%    Version: 22-Oct-2024
%    Comments by ChatGPT4, 13-May-2025
%  See also:
%    - Section 12.2 of *EIFEM_largeROTfinal.pdf*
%    - Bmat_weights_EIFEbub.m (non-corotational version)
%    - EIFE_CONSISTENCY.m (for validation)
%
%--------------------------------------------------------------------------




% Modification of Bmat_weights_EIFEbub.m for corotational formulation
% JAHO, 22-oct-2024, Balmes 185, Barcelona
% See  /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/02_PURE_ROTATION.mlx
if nargin == 0
    load('tmp3.mat')
end

if DATA.UNIFORM_SCALING_REFERENCE_ELEMENT == 1
    indCHOSEN = TRANSF_COORD.IndexParentDomain ;
    
    % Gradient matrices 
   %  \BmatIstE{e} \defeq  \dfrac{1}{\lambdaLENe{e}} \BmatIrefSTe{e}, \hspace{0.5cm} e = 1,2 \ldots \nelemC
    BmatIstE = EIFEoper_all(indCHOSEN).INTforces.BmatIrefST/TRANSF_COORD.lambdaLEN ;
    
    % CECM WEIGHTS
    %  \WecmE{e} = {\lambdaLENe{e}}^{\nSD}  \WecmREFe{e}, \hspace{0.5cm} e = 1,2 \ldots \nelemC
    WecmE = TRANSF_COORD.detJe*EIFEoper_all(indCHOSEN).INTforces.weights ;  % CECM weights
  
    
    
    % Linear stiffness matrix
    % ************************
    %  {\KcLINlocE{e}} = {\lambdaLENe{e}}^{\nSD-2}   {\KcLINrefE{e}}
    nSD = length(TRANSF_COORD.CentCglo) ;
    lambda_n2 = TRANSF_COORD.lambdaLEN^(nSD-2) ;
    KcLINlocE= lambda_n2*EIFEoper_all(indCHOSEN).Kcoarse  ;

 
    
    %     if ~isempty(EIFEoper_all(indCHOSEN).INFO.DOFsBUB)
    %         % There are bubble modes
    %         % --------------------------------------------
    %         nmodesALL = size(EIFEoper_all(indCHOSEN).INTforces.BmatRED,2); % This is the total number of modes
    %         ndofsBUB = length(EIFEoper_all(indCHOSEN).INFO.DOFsBUB) ;  % This is the number of bubble modes
    %         indBOUND = 1:(nmodesALL-ndofsBUB) ;   % These are the indexes of the strain modes
    %         indBUB = (nmodesALL-ndofsBUB+1):nmodesALL  ; % These are the indexes of the bubble modes
    %         %
    
    
    %  BmatB = (EIFEoper_all(indCHOSEN).INTforces.BmatRED(:,indBOUND)*EIFEoper_all(indCHOSEN).OPER.HdefINV_PsiDEFfT*Vrot)/TRANSF_COORD.SCALEFACTOR;
    %  BmatBUB = EIFEoper_all(indCHOSEN).INTforces.BmatRED(:,indBUB)/TRANSF_COORD.SCALEFACTOR ;
    
    % Extended B-matrix
    % --------------------------
    %  ndim = size(TRANSF_COORD.ROTATION,2) ;
    %  Bmat = ExtendedBmatrixBUB(BmatB,BmatBUB,ndim,ndofsBUB) ;
    
    % Elastic coarse-scale stiffness
    
    
    %     else
    %         Bmat =  (EIFEoper_all(indCHOSEN).INTforces.BmatRED*EIFEoper_all(indCHOSEN).OPER.HdefINV_PsiDEFfT*Vrot)/TRANSF_COORD.SCALEFACTOR;
    %         % Elastic coarse-scale stiffness
    %         ndofsBUB = 0 ;
    %         KcELAS = ElasticCoarseStiff_1elem(EIFEoper_all,indCHOSEN,TRANSF_COORD,ndofsBUB) ;
    %     end
else
    error('Option not implemented yet')
    % Future implementations (11-March-2023) should consider the
    % preliminary version developed
    %  /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/01_HOMOG/02_POST_PROCESS/NonVectorElastCode/ComputeKeMatrix_DilRot_multi.m
end

