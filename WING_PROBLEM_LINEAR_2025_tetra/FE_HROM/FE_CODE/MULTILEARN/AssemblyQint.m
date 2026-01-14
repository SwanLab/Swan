function  [J]=  AssemblyQint(BasisRdef,NODESbound,reactDOMrbGLO,betaBC,COL,COLloc,...
    THETAfaces,GAMMAentities,ndim,alphaBC,r,BasisRrb,BasisUdef,BasisUrb,COORref)

%dbstop('5')
if nargin == 0
    load('tmp.mat')
    %GAMMAentities = GAMMAfaces ;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(BasisRdef) >1
    error('Option not implemented for several types of RVE')
end
nRB= size(BasisRrb,2) ;
nDOMrve = cellfun(@length,COL) ;   % Number of domains
[nDOF nMODES ]= cellfun(@size,BasisRdef) ;  % Number of reaction modes
nMODESreac = nMODES + nRB ;
nrows  = sum((nDOMrve+1).*nRB) ;
ncols = sum(nDOMrve.*nMODESreac) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
J = sparse(nrows,ncols) ;
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
f1 = NODESfaces{1} ;
COORf1  = COORref(f1,:) ;
% REFERENCE POINT
COORpoint = sum(COORf1,1)/size(COORf1,1) ;
COORreference = bsxfun(@plus,COORf1',-COORpoint')';

Q = ConstructBasisRigidBody(COORreference) ;


f1 = small2large(NODESfaces{1},ndim) ;
f2 = small2large(NODESfaces{3},ndim) ;






iROWS = 0 ;
iCOLS = 0 ;
nMODES_all = zeros(1,nDOM) ;
for itype = 1:length(BasisRdef)
    nMODES_all(COL{itype}) = size(BasisRdef{itype},2) ;
end
BasisR = [BasisRrb,BasisRdef{1}] ;
QR_f1 = Q'*BasisR(f1,:) ;
QR_f2 = Q'*BasisR(f2,:) ;



for iface = 1:nDOM+1
    if iface == 1
        
        VECT=  [(1-ALPHA(1))*QR_f1  ] ;
        indROWS = 1:size(VECT,1) ;
        indCOLS = 1:size(VECT,2) ;
        J(indROWS,indCOLS) = VECT ;
        iROWS = length(indROWS) ; iCOLS = 0 ;
    elseif iface == (nDOM+1)
        VECT=  [(1-ALPHA(2))*QR_f2  ] ;
        indROWS = iROWS + (1:size(VECT,1)) ;
        indCOLS = iCOLS + (1:size(VECT,2))  ;
        J(indROWS,indCOLS) =  VECT ;
    else
       
        VECT = [QR_f2, QR_f1] ;
        indROWS = iROWS + (1:size(VECT,1)) ;
        indCOLS = iCOLS + (1:size(VECT,2))  ;
        J(indROWS,indCOLS) = VECT;
        
        iROWS = iROWS + length(indROWS) ;
         iCOLS = iCOLS + size(QR_f2,2) ; 
        
    end
end



end






