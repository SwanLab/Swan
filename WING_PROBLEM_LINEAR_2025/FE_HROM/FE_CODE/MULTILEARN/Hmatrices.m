function Hqr_d = Hmatrices(BasisU,BasisR,f1,f2,nDOM,ALPHA)
Cov_f1f1 = sparse(BasisU(f1,:)'*BasisR(f1,:)) ;
Cov_f2f2 = sparse(BasisU(f2,:)'*BasisR(f2,:)) ;
Cov_f= sparse(Cov_f1f1+Cov_f2f2) ;

if nDOM==1
    
    if ALPHA(1) ==1 & ALPHA(2) ==0
        Hqr_d =Cov_f1f1 ;
    elseif ALPHA(2) ==1 & ALPHA(1) ==0
        Hqr_d =Cov_f2f2 ;
    elseif  ALPHA(2) ==1 & ALPHA(1) ==1
        Hqr_d =Cov_f;
    else
        error('Ill-posed problem')
    end
    
else
    
    Hqr_d = cell(nDOM,1) ;
    Hqr_d(:) = {Cov_f} ;
    Hqr_d{1} = ALPHA(1)*Cov_f1f1 +Cov_f2f2   ;
    Hqr_d{end} = ALPHA(2)*Cov_f2f2 +Cov_f1f1   ;
    Hqr_d = blkdiag(Hqr_d{:}) ;
end