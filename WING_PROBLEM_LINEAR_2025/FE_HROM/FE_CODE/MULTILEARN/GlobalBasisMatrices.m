function [BasisUrbGLO,BasisUdefGLO,BasisRdefGLO] = GlobalBasisMatrices(BasisUrb,nDOM,BasisUdef,BasisRdef,COLUMNS_RVE)

if nargin == 0
    load('tmp.mat')
end

BasisUrbGLO = repmat(sparse(BasisUrb),1,nDOM) ;
BasisUrbGLO = mat2cell(BasisUrbGLO,size(BasisUrb,1),repmat(size(BasisUrb,2),nDOM,1)) ;
BasisUrbGLO = blkdiag(BasisUrbGLO{:}) ;
% Deformation basis matrices
BasisUdefGLO = DiagonalGlobalMatrixRVEs(nDOM,BasisUdef,COLUMNS_RVE) ;
% Reaction basis matrices
BasisRdefGLO = DiagonalGlobalMatrixRVEs(nDOM,BasisRdef,COLUMNS_RVE) ;