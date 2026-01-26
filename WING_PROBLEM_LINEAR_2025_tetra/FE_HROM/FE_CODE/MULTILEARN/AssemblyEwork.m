function  [E]=  AssemblyEwork(BasisRdef,NODESbound,reactDOMrbGLO,betaBC,COL,COLloc,...
    THETAfaces,GAMMAentities,ndim,alphaBC,r,BasisRrb,BasisUdef,BasisUrb)

%dbstop('5')
if nargin == 0
    load('tmp.mat')
    %GAMMAentities = GAMMAfaces ;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nRB= size(BasisRrb,2) ;
nDOMrve = cellfun(@length,COL) ;   % Number of domains
[nDOF nMODES ]= cellfun(@size,BasisRdef) ;  % Number of reaction modes
nMODESreac = nMODES + nRB ;
nrowsREAC  = sum(nDOMrve.*nMODESreac) ;
[nDOF nMODES ]= cellfun(@size,BasisUdef) ;  % Number of reaction modes
nMODESdef = nMODES + nRB ;
ncolsDEF  = sum(nDOMrve.*nMODESdef) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E = sparse(nrowsREAC,ncolsDEF) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nDOM = size(betaBC,1) ;
nMODES_all = zeros(1,nDOM) ;
for itype = 1:length(BasisRdef)
    nMODES_all(COL{itype}) = size(BasisRdef{itype},2) + nRB ;
end
INDrig = [] ;
iacum = 0 ;
 
ALPHA = zeros(1,2) ; 
if alphaBC(1,1) ==1 
ALPHA(1) = 1; 
end
if alphaBC(end,3) ==1 
ALPHA(2) = 1 ; 
end
 
NODESfaces = NODESbound.PLANE ; 
f1 = small2large(NODESfaces{1},ndim) ; 
f2 = small2large(NODESfaces{3},ndim) ; 

iROWS = 0 ;
iCOLS = 0 ;
nMODES_all = zeros(1,nDOM) ;
for itype = 1:length(BasisRdef)
    nMODES_all(COL{itype}) = size(BasisRdef{itype},2) ;
end
for idom = 1:nDOM
    itype_i = COLloc(idom) ;  % Type of RVE
    BasisRi = [BasisRrb,BasisRdef{itype_i}] ;     
    BasisUi = [BasisUrb,BasisUdef{itype_i}] ;   
    Cov_f1f1 = BasisRi(f1,:)'*BasisUi(f1,:) ;
    Cov_f2f2 = BasisRi(f2,:)'*BasisUi(f2,:) ;
    if idom == 1
        if nDOM >1            
            itype_ip1 = COLloc(idom+1) ;
            BasisRip1 = [BasisRrb,BasisRdef{itype_ip1}] ;
            BasisUip1 = [BasisUrb,BasisUdef{itype_ip1}] ;
            Cov_f2f1 = BasisRi(f2,:)'* BasisUip1(f1,:) ;
            VECT = [(1-ALPHA(1))*Cov_f1f1 +  Cov_f2f2, -Cov_f2f1] ;
            indROWS = 1:size(VECT,1) ;
            indCOLS = 1:size(VECT,2) ;
            E(indROWS,indCOLS) = VECT ;
            iROWS = length(indROWS) ; iCOLS = 0 ; 
        else
            VECT=  [(1-ALPHA(1))*Cov_f1f1 + (1-ALPHA(end))*Cov_f2f2] ; 
            indROWS = 1:size(VECT,1) ;
            indCOLS = 1:size(VECT,2) ;
            E(indROWS,indCOLS) = VECT ; 
         end
    elseif idom == nDOM
        itype_im1 = COLloc(idom-1) ; % Type of contiguous RVE
        BasisRim1 = [BasisRrb,BasisRdef{itype_im1}] ;
        BasisUim1 = [BasisUrb,BasisUdef{itype_im1}] ;
        Cov_f1f2 = BasisRi(f1,:)'* BasisUim1(f2,:) ;
         VECT = [-Cov_f1f2, Cov_f1f1+ (1-ALPHA(end))*Cov_f2f2] ; ;
        indROWS = iROWS + (1:size(VECT,1)) ;
        indCOLS = iCOLS + (1:size(VECT,2))  ;
        E(indROWS,indCOLS) =  VECT ;       
    else
        itype_im1 = COLloc(idom-1) ; % Type of contiguous RVE
        itype_ip1 = COLloc(idom+1) ; % Type of contiguous RVE
        BasisRim1 = [BasisRrb,BasisRdef{itype_im1}] ;
        BasisUim1 = [BasisUrb,BasisUdef{itype_im1}] ;
        BasisRip1 = [BasisRrb,BasisRdef{itype_ip1}] ;
        BasisUip1 = [BasisUrb,BasisUdef{itype_ip1}] ;

        Cov_f2f1 = BasisRi(f2,:)'*BasisUip1(f1,:) ;
        Cov_f1f2 = BasisRi(f1,:)'*BasisUim1(f2,:) ;%
        VECT = [+Cov_f1f2, -Cov_f1f1 + Cov_f2f2  , -Cov_f2f1] ;
        indROWS = iROWS + (1:size(VECT,1)) ;
        indCOLS = iCOLS + (1:size(VECT,2))  ;
        E(indROWS,indCOLS) = VECT;
       
        iROWS = iROWS + length(indROWS) ;
        iCOLS = sum(nMODES_all(1:idom-1)) ;
        
    end
end



end






