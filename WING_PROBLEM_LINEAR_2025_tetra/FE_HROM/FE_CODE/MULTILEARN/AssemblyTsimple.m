function  [T]=  AssemblyTsimple(BasisRdef,NODESbound,reactDOMrbGLO,betaBC,COL,COLloc,...
    THETAfaces,GAMMAentities,ndim,alphaBC,r,BasisRrb,BasisUdef,BasisUrb,COORref,uBAR)

%dbstop('5')
if nargin == 0
    load('tmp.mat')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(BasisRdef) >1
    error('Option not implemented for several types of RVE')
end
nRB= size(BasisRrb,2) ;
nDOMrve = cellfun(@length,COL) ;   % Number of domains
[nDOF nMODES ]= cellfun(@size,BasisRdef) ;  % Number of reaction modes
nMODESreac = nMODES + nRB ;
ncols  = sum((nDOMrve).*nMODESreac) ;
[nDOF nMODES ]= cellfun(@size,BasisUdef) ;  % Number of reaction modes
nMODESdef = nMODES + nRB ;

nrows = sum((nDOMrve+1).*nMODESdef) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T = sparse(nrows,ncols) ;
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

%
 f1 = small2large(NODESfaces{1},ndim) ;
f2 = small2large(NODESfaces{3},ndim) ;






iROWS = 0 ;
iCOLS = 0 ;
nMODES_all = zeros(1,nDOM) ;
for itype = 1:length(BasisRdef)
    nMODES_all(COL{itype}) = size(BasisRdef{itype},2) ;
end
BasisR = [BasisRrb,BasisRdef{1}] ;
BasisU = [BasisUrb,BasisUdef{1}] ;

Cov_f1f1 = BasisU(f1,:)'*BasisR(f1,:) ; 
Cov_f1f2 = BasisU(f1,:)'*BasisR(f2,:) ; 
Cov_f2f2 = BasisU(f2,:)'*BasisR(f2,:) ; 




for iface = 1:nDOM+1
    if iface == 1
        
        VECT=  [ (1-ALPHA(1))*Cov_f1f1  ] ;
        indROWS = 1:size(VECT,1) ;
        indCOLS = 1:size(VECT,2) ;
        T(indROWS,indCOLS) = VECT ;
        iROWS = length(indROWS) ; iCOLS = 0 ;
    elseif iface == (nDOM+1)
        VECT=  [(1-ALPHA(2))*Cov_f2f2  ] ;
        indROWS = iROWS + (1:size(VECT,1)) ;
        indCOLS = iCOLS + (1:size(VECT,2))  ;
        T(indROWS,indCOLS) =  VECT ;
    else
       
        VECT = [+Cov_f1f2, Cov_f1f1] ;
        indROWS = iROWS + (1:size(VECT,1)) ;
        indCOLS = iCOLS + (1:size(VECT,2))  ;
        T(indROWS,indCOLS) = VECT;
        
        iROWS = iROWS + length(indROWS) ;
         iCOLS = iCOLS + size(Cov_f1f2,2) ; 
        
    end
end



end






