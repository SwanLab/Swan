function  [ Prr Prd Pdd bCOMPr bCOMPd ]=  AssemblyPcomp(BasisUrb,BasisUdef,f1,f2,nDOM,uBAR)

%dbstop('4')
BasisUdom = [BasisUrb BasisUdef] ;
f = [f1 ; f2] ;
Cov_f = BasisUdom(f,:)'*BasisUdom(f,:) ;
Cov_f2f1 = BasisUdom(f2,:)'*BasisUdom(f1,:) ;
Cov_f1f2 = BasisUdom(f1,:)'*BasisUdom(f2,:) ;
Cov_f1f1 = BasisUdom(f1,:)'*BasisUdom(f1,:) ;
nMODES = size(BasisUdom,2) ;
nrows = nMODES*nDOM ;
P = sparse(nrows,nrows) ;
nRB = size(BasisUrb,2) ;
nDEF = nMODES-nRB ;
INDrig = [] ;

for idom = 1:nDOM
    if idom == 1
        if nDOM >1
            indROWS = 1:nMODES ;
            indCOLS = 1:2*nMODES ;
            P(indROWS,indCOLS) = [Cov_f, -Cov_f2f1] ;
        else
            indROWS = 1:nMODES ;
            indCOLS = 1:nMODES ;
            P(indROWS,indCOLS) = [Cov_f1f1] ;
        end
    elseif idom == nDOM
        indROWS = ((nDOM-1)*nMODES+1):nDOM*nMODES ;
        indCOLS = ((nDOM-2)*nMODES+1):nDOM*nMODES ;
        P(indROWS,indCOLS) = [-Cov_f1f2, Cov_f1f1] ;
    else
        indROWS = ((idom-1)*nMODES+1):idom*nMODES ;
        indCOLS = ((idom-2)*nMODES+1):(idom+1)*nMODES ;
        P(indROWS,indCOLS) = [-Cov_f1f2, Cov_f, -Cov_f2f1] ;
    end
    INDrig =[INDrig;[(idom-1)*nMODES+(1:nRB)]'] ;
end

INDdef = (1:nrows)';
INDdef(INDrig) =  [] ;
Prr = P(INDrig,INDrig) ;
Prd = P(INDrig,INDdef) ;
Pdd = P(INDdef,INDdef) ;
bCOMPr = zeros(size(INDrig)) ;
indROWS = 1:nRB ;
bCOMPr(indROWS) = BasisUrb(f1,:)'*uBAR ;
bCOMPd = zeros(size(INDdef)) ;

indROWS = 1:nDEF ;
bCOMPd(indROWS) = BasisUdef(f1,:)'*uBAR ;




