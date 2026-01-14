function  [T, c, s]=  AssemblyTreac(BasisRdef,f1,f2,nDOM,reactDOMrb)

if nargin == 0
    load('tmp2.mat')
end


f = [f1 ; f2] ;
Cov_f = BasisRdef(f,:)'*BasisRdef(f,:) ;
Cov_f2f1 = BasisRdef(f2,:)'*BasisRdef(f1,:) ;
Cov_f1f2 = BasisRdef(f1,:)'*BasisRdef(f2,:) ;
Cov_f2f2 = BasisRdef(f2,:)'*BasisRdef(f2,:) ;
nMODES = size(BasisRdef,2) ;
nrows = nMODES*nDOM ;
T = sparse(nrows,nrows) ;
c = zeros(nrows,1);
s = 0 ;
for idom = 1:nDOM
    if idom == 1
        if nDOM >1
            indROWS = 1:nMODES ;
            indCOLS = 1:2*nMODES ;
            T(indROWS,indCOLS) = [Cov_f2f2, Cov_f2f1] ;
            c(indROWS) = (reactDOMrb{idom}(f2,:) + reactDOMrb{idom+1}(f1,:))'*BasisRdef(f2,:) ;
            s = s + norm(reactDOMrb{idom}(f2,:) + reactDOMrb{idom+1}(f1,:))^2 ;
        else
            indROWS = 1 ;
            indCOLS = 1;
            T(indROWS,indCOLS) = [Cov_f2f2] ;
            c(indROWS) = (reactDOMrb{idom}(f2,:) )'*BasisRdef(f2,:) ;
            s = s + norm(reactDOMrb{idom}(f2,:))^2 ;
        end
    elseif idom == nDOM
        indROWS = ((nDOM-1)*nMODES+1):nDOM*nMODES ;
        indCOLS = ((nDOM-2)*nMODES+1):nDOM*nMODES ;
        T(indROWS,indCOLS) = [Cov_f1f2, Cov_f] ;
        c(indROWS) = (reactDOMrb{idom-1}(f2,:) + reactDOMrb{idom}(f1,:))'*BasisRdef(f1,:) + ...
            reactDOMrb{idom}(f2,:)'*BasisRdef(f2,:);
        s = s + norm(reactDOMrb{idom}(f2,:))^2 ;
    else
        indROWS = ((idom-1)*nMODES+1):idom*nMODES ;
        indCOLS = ((idom-2)*nMODES+1):(idom+1)*nMODES ;
        T(indROWS,indCOLS) = [Cov_f1f2, Cov_f, Cov_f2f1] ;
        c(indROWS) = (reactDOMrb{idom-1}(f2,:) + reactDOMrb{idom}(f1,:))'*BasisRdef(f1,:) + ...
            (reactDOMrb{idom}(f2,:) + reactDOMrb{idom+1}(f1,:))'*BasisRdef(f2,:)   ;
        s = s + norm(reactDOMrb{idom}(f2,:) + reactDOMrb{idom+1}(f1,:))^2 ;
        
    end
end









