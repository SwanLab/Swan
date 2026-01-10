function Tint =  Tmatrix(BasisR,f1,f2,nDOM,BasisINT,ISTD,DATAIN)
Covf2 = (BasisINT'*BasisR(f2,:))';
Covf1 = (BasisINT'*BasisR(f1,:))' ;

norm_columnsF1 =  sqrt(sum(Covf1.^2,1)) ; 
norm_columnsF2 =  sqrt(sum(Covf2.^2,1)) ; 

norm_rowsF1 =  sqrt(sum(Covf1.^2,2)) ; 
norm_rowsF2 =  sqrt(sum(Covf2.^2,2)) ; 

if ISTD == 1
    DATAIN = DefaultField(DATAIN,'TOLERANCE_BasisINT_times_BasisRdef_SVD',1e-4) ;
    TOL =  DATAIN.TOLERANCE_BasisINT_times_BasisRdef_SVD ;
    [dummy sF2] = SVDT(Covf2,0) ;
    [dummy sF1]  = SVDT(Covf1,0) ;
    
    
    if  (sF2(end)) < TOL
        disp(['BasisINT*BasisRdef_f2 ='])
        disp(num2str(Covf2))
        sF2
        error(['Rank deficient Tint matrix. Enrich the training set of increase the number of disp. modes'])
    end
    if  (sF1(end)) < TOL
        disp(['BasisINT*BasisRdef_f1 ='])
        disp(num2str(Covf1))
        sF1
        error(['Rank deficient Tint matrix. Enrich the training set of increase the number of disp. modes'])
    end
    
end

%% "Displaced" Sum of diagonal matrices
if nDOM>1
    DIAG1 = cell(nDOM-1,1) ;
    DIAG1(:) = {sparse(Covf2)} ;
    DIAG1 = blkdiag(DIAG1{:}) ;
    DIAG2 = cell(nDOM-1,1) ;
    DIAG2(:) = {sparse(Covf1)} ;
    DIAG2 = blkdiag(DIAG2{:}) ;
    
    ncols = (nDOM-1)*size(BasisINT,2) ;
    MATZEROS = sparse(size(BasisR,2),ncols) ;
    
    Tint = [DIAG1;MATZEROS] + [MATZEROS;DIAG2] ;
    
else
    %  error('Number of domains must be greater than 1')
    Tint = [] ;
end
Tint = Tint' ;