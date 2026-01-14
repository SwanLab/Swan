function [wSTs Bst, BstW, wST,  DATAIN] = ...
    ComputeBmatricesROM(DATAROM_glo,DATA_REFMESH_glo,DATAIN,ndim,MESH1D)
if nargin == 0
    load('tmp0.mat')
end

% Computing Coarse-scale B-matrices
% -----------------------------------
nelem = length(DATAROM_glo) ;
Belem = cell(nelem,1) ;
wSTs = cell(nelem,1) ;
wST = cell(nelem,1) ;

nstrain = DATAROM_glo{1}.HROMVAR.nstrain ; 

for itype = 1:length(DATAROM_glo)
    % This loop is actually not needed ---there is just one type of SLICE
    BCdom = DATAROM_glo{itype}.HROMVAR.BCdom ;
    ELEMS = find(MESH1D.MaterialType == itype) ;
    Belem(ELEMS,1) = {BCdom} ;
    wSTs(ELEMS,1) = {DATAROM_glo{itype}.HROMVAR.WdomRED} ; 
    Wrep = repmat(DATAROM_glo{itype}.HROMVAR.WdomRED',nstrain,1) ;
    %Wrep = sparse(Wrep(:)) ;
    %Wdiag = diag(Wrep) ; 
    wST(ELEMS,1) = {Wrep(:)} ; 
    
end

Belem = cell2mat(Belem) ; 
wSTs = cell2mat(wSTs) ; 
wST  = cell2mat(wST) ; 
ngaus = length(DATAROM_glo{1}.HROMVAR.setPoints); 
CN = MESH1D.CN ; 
nnode = length(unique(CN)) ; 

nnodeE = 2 ; 
nelem = size(CN,1) ; 

Bst = AssemblyBGlobal(Belem,nstrain,nelem,nnodeE,ndim,ngaus,CN,nnode) ;

wDIAG = diag(sparse(wST)) ; 
BstW = wDIAG*Bst ; 

% Matrix multiplied by weights 
 
