function  [ Prr Prd Pdd bCOMPr bCOMPd ]=  AssemblyPcompMULTI(BasisUrb,BasisUdef,f1,f2,nDOM,uBAR,ALPHA,COL,COLloc)
%dbstop('3')
if nargin == 0
    load('tmp.mat')
end
f = [f1 ; f2] ;
nRB = size(BasisUrb,2) ; % Number of rigid body modes
nDOMrve = cellfun(@length,COL) ;   % Number of domains
[nDOF nDEFrve ]= cellfun(@size,BasisUdef) ;  % Number of deformation modes
nMODESrve = nDEFrve + nRB ; % Total number of modes
nrows  = sum(nDOMrve.*nMODESrve) ;
P = sparse(nrows,nrows) ;
INDrig = [] ;
iacum = 0 ;

nMODES_all = zeros(1,nDOM) ;
for itype = 1:length(BasisUdef)
    nMODES_all(COL{itype}) = size(BasisUdef{itype},2) + nRB ;
end
iROWS =0 ; iCOLS = 0 ;
for idom = 1:nDOM
    itype_i = COLloc(idom) ;  % Type of RVE
    BasisUdom_i = [BasisUrb BasisUdef{itype_i}] ;  % Basis matrix of displacements for this type of RVE
    Cov_f1f1 = BasisUdom_i(f1,:)'*BasisUdom_i(f1,:) ;
    Cov_f2f2 = BasisUdom_i(f2,:)'*BasisUdom_i(f2,:) ;
    Cov_f = BasisUdom_i(f,:)'*BasisUdom_i(f,:) ;
    if idom == 1        
        if nDOM >1
            itype_ip1 = COLloc(idom+1) ; % Type of contiguous RVE
            BasisUdom_ip1 = [BasisUrb BasisUdef{itype_ip1}] ;
            Cov_f2f1 = BasisUdom_i(f2,:)'*BasisUdom_ip1(f1,:) ;
            indROWS = 1:size(BasisUdom_i,2) ;
            indCOLS = 1:(size(BasisUdom_i,2)+ size(BasisUdom_ip1,2));
            P(indROWS,indCOLS) = [ALPHA(1)*Cov_f1f1+Cov_f2f2, -Cov_f2f1] ;
            iROWS = length(indROWS) ; iCOLS = 0 ; 
        else
            indROWS = 1:size(BasisUdom_i,2)  ;
            indCOLS = 1:size(BasisUdom_i,2) ;
            P(indROWS,indCOLS) = [ALPHA(1)*Cov_f1f1 + ALPHA(2)*Cov_f2f2] ;
        end
    elseif idom == nDOM
        itype_im1 = COLloc(idom-1) ; % Type of contiguous RVE
        BasisUdom_im1 = [BasisUrb BasisUdef{itype_im1}] ;
        Cov_f1f2= BasisUdom_i(f1,:)'*BasisUdom_im1(f2,:) ;        
        indROWS = iROWS + (1:size(Cov_f1f2,1)) ;
        indCOLS = iCOLS + (1:(size(Cov_f1f2,2)+size(Cov_f1f1,2) ))  ;
        P(indROWS,indCOLS) = [-Cov_f1f2, Cov_f1f1 + ALPHA(2)*Cov_f2f2] ;
    else
        itype_im1 = COLloc(idom-1) ; % Type of contiguous RVE
        BasisUdom_im1 = [BasisUrb BasisUdef{itype_im1}] ;
        itype_ip1 = COLloc(idom+1) ; % Type of contiguous RVE
        BasisUdom_ip1 = [BasisUrb BasisUdef{itype_ip1}] ;
        Cov_f2f1 = BasisUdom_i(f2,:)'*BasisUdom_ip1(f1,:) ;
        Cov_f1f2= BasisUdom_i(f1,:)'*BasisUdom_im1(f2,:) ;
        
        indROWS = iROWS + (1:size(Cov_f1f2,1)) ;
        indCOLS = iCOLS + (1:(size(Cov_f1f2,2)+size(Cov_f,2) )+size(Cov_f2f1,2))  ;
        P(indROWS,indCOLS) = [-Cov_f1f2, Cov_f, -Cov_f2f1] ;
        iROWS = iROWS + length(indROWS) ; 
        iCOLS = sum(nMODES_all(1:idom-1)) ;
    end
    %¡  INDrig =[INDrig;[(idom-1)*nMODES+(1:nRB)]'] ;
    INDrig =[INDrig;[iacum+(1:nRB)]'] ;
    iacum =  sum(nMODES_all(1:idom)) ;
end

INDdef = (1:nrows)';
INDdef(INDrig) =  [] ;
Prr = P(INDrig,INDrig) ;
Prd = P(INDrig,INDdef) ;
Pdd = P(INDdef,INDdef) ;

% Rigid body vector
bCOMPr = zeros(size(INDrig)) ;
indROWS = 1:nRB ;
bCOMPr(indROWS) = ALPHA(1)*BasisUrb(f1,:)'*uBAR{1} ;
indROWS = (length(bCOMPr)-nRB+1):length(bCOMPr) ;
bCOMPr(indROWS) = ALPHA(2)*BasisUrb(f2,:)'*uBAR{2} ;

% Deformation vector
bCOMPd = zeros(size(INDdef)) ;
itype_i = COLloc(idom) ;  % Type of RVE
nDEF = size(BasisUdef{itype_i},2) ; 
indROWS = 1:nDEF ;
bCOMPd(indROWS) = ALPHA(1)*BasisUdef{itype_i}(f1,:)'*uBAR{1} ;
itype_i = COLloc(end) ;  % Type of RVE
nDEF = size(BasisUdef{itype_i},2) ; 
indROWS = (length(bCOMPd)-nDEF+1):length(bCOMPd) ;
bCOMPd(indROWS) = ALPHA(2)*BasisUdef{itype_i}(f2,:)'*uBAR{2} ;




