function [Hqr KdomREDglo] = RedStiffnessMatrix_Hmatrix(BasisUdef,nDOM,COLUMNS_RVE,KdomRED,BasisRdef,...
    DATAINM,BasisUrb,BasisRrb)

if nargin == 0
    load('tmp2.mat')
end


if DATAINM.MinimizationBoundaryWork == 0 | DATAINM.MinimizationBoundaryWork == 2
    COV = cell(size(BasisUdef)) ;
    for itype = 1:length(BasisUdef)
        COV{itype} =BasisUdef{itype}'*BasisRdef{itype} ;
    end
    Hqr = DiagonalGlobalMatrixRVEs(nDOM,COV,COLUMNS_RVE) ;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % E) Stiffness matrix
    % -----------------------------
    KdomREDglo = DiagonalGlobalMatrixRVEs(nDOM,KdomRED,COLUMNS_RVE) ;
elseif DATAINM.MinimizationBoundaryWork == 1
    
    COV = cell(size(BasisUdef)) ;
    for itype = 1:length(BasisUdef)
        BasisU = [BasisUrb  BasisUdef{itype}] ;
        BasisR = [BasisRrb  BasisRdef{itype}] ;
        COV{itype} =BasisU'*BasisR ;
    end
    Hqr = DiagonalGlobalMatrixRVEs(nDOM,COV,COLUMNS_RVE) ;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % E) Stiffness matrix
    % -----------------------------
    COV = cell(size(BasisUdef)) ;
    for itype = 1:length(BasisUdef)
        nrows = size(BasisUrb,2) + size(BasisUdef{itype},2) ;
        COV{itype} = sparse(nrows,nrows) ;
        irows = (size(BasisUrb,2) +1):nrows ;
        COV{itype}(irows,irows) = KdomRED{itype} ;
    end
    KdomREDglo = DiagonalGlobalMatrixRVEs(nDOM,COV,COLUMNS_RVE) ;
    
else
    error('option not implemented')
end