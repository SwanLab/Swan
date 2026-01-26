function BasisUdefGLO = DiagonalGlobalMatrixRVEs(nDOM,BasisUdef,COLUMNS_RVE)


BasisUdefGLO = cell(1,nDOM) ; 
for itype = 1:length(COLUMNS_RVE)
    COL = COLUMNS_RVE{itype} ;
    MATLOC = repmat(sparse(BasisUdef{itype}),1,length(COL)) ;  
    MATLOC = mat2cell(MATLOC,size(BasisUdef{itype},1),repmat(size(BasisUdef{itype},2),length(COL),1)) ;
    BasisUdefGLO(COL) = MATLOC ;  
end
%BasisUdefGLO = mat2cell(BasisUdefGLO,size(BasisUdef,1),repmat(size(BasisUdef,2),nDOM,1)) ;
BasisUdefGLO = blkdiag(BasisUdefGLO{:}) ;