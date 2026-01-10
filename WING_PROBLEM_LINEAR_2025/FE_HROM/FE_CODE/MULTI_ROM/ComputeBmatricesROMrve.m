function [wSTs Bst, BstW, wST,  DATAIN,DOFsKEEP,ASSEMBLY_INFO] = ...
    ComputeBmatricesROMrve(DATAROM_glo,DATA_REFMESH_glo,DATAIN,ndimALL,MESH2D)
if nargin == 0
    load('tmp2.mat')
end


%% We have that ndim(1) = ndim(3), and ndim(2) = 4
%  The idea is to make ndimMAX = max(ndim), and proceed as in the normal
%  case, but filling Kskel with zeros (so that it is 4 ndimMAX x 4 ndimMAX)  
% ------------------------------------------------------------------
ndim = max(ndimALL) ;
[DOFSeliminate ]= DOFSeliminateNEQDIM(ndimALL,MESH2D) ; 


% % ---------------------------------------------
% for itype = 1:length(DATAROM)
%     Ke = DATAROM{itype}.Kskel ;
%     Ke = mat2cell(Ke,ndim,ndim); 
%     for iface = 1:size(Ke,1)
%         for jface =1:size(Ke,2) 
%             [nrows,ncols] = size(Ke{iface,jface}) ;  
%             KeLOCAL = zeros(ndimMAX,ndimMAX) ;
%             KeLOCAL(1:nrows,1:ncols) = Ke{iface,jface} ; 
%             Ke{iface,jface} = KeLOCAL ; 
%         end
%     end
%      DATAROM{itype}.Kskel = cell2mat(Ke);    
% end

% Computing Coarse-scale B-matrices
% -----------------------------------
nelem = length(DATAROM_glo) ;
Belem = cell(nelem,1) ;
wSTs = cell(nelem,1) ;
wST = cell(nelem,1) ;

nstrain = DATAROM_glo{1}.HROMVAR.nstrain ; 

for itype = 1:length(DATAROM_glo)
    % This loop is actually not needed ---there is just one type of SLICE
    % ----------------------------------------------------------------------
    BCdom = DATAROM_glo{itype}.HROMVAR.BCdom ; %  This is the B matrix 
    % of the coarse scale model. It features as many rows as integrations
    % points times number of strain components, and as many columns as
    % coarse-scale DOFs.  Notice that 
    % BCdom = [BCdom_1, BCdom_2, BCdom_3, BCdom_4]. Here, the index refers
    % to the face of the RVE. 
    nfaces= length(ndimALL) ; 
  %  BCdom_extended = zeros(size(BCdom,1),ndim*nfaces) ; % Extended B matrix 
        BCdom_extended = zeros(size(BCdom,1),ndim*nfaces) ; % Extended B matrix ....  

    % All faces are assumed to have the same number of DOFs
    iiniOLD = 1; iiniNEW = 1;  
    for iface = 1:nfaces
        ifinOLD = iiniOLD+ndimALL(iface)-1 ; 
        ifinNEW =  iiniNEW+ndim-1 ;  ifinNEWOLD = iiniNEW + ndimALL(iface)-1 ; 
        BCdom_extended(:,iiniNEW:ifinNEWOLD) =  BCdom(:,iiniOLD:ifinOLD) ;
        iiniOLD = ifinOLD +1;  iiniNEW = ifinNEW +1 ; 
    end
    
    ELEMS = find(MESH2D.MaterialType == itype) ;
    Belem(ELEMS,1) = {BCdom_extended} ;
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
CN = MESH2D.CN ; 
nnode = length(unique(CN)) ; 

nnodeE = size(CN,2) ; 
nelem = size(CN,1) ; 

Bst = AssemblyBGlobal(Belem,nstrain,nelem,nnodeE,ndim,ngaus,CN,nnode) ;

wDIAG = diag(sparse(wST)) ; 

DATA_INCLUDE_WEIGHTS = 0 ;  % 25-June-2019 -.... For the sake of efficiency !!! 
if DATA_INCLUDE_WEIGHTS == 1
    error('This option is memory demanding ... Abandoned on 25-June-2019')
    BstW = wDIAG*Bst ;
else
    BstW  = [] ;
end


% % Now we remove the indices DOFSeliminate 
% % -----------------------------------------
DOFSeliminate = sort(DOFSeliminate) ; 
DOFsKEEP = setdiff(1:size(Bst,2),DOFSeliminate) ; 
Bst = Bst(:,DOFsKEEP) ; 
if ~isempty(BstW)
BstW = BstW(:,DOFsKEEP) ; 
end


% For the case of standard vectorize assembly (no  B*C*B method)
if  DATAIN.AssemblyMethodK_BCB == 0
    error('Method not implemented')
    ASSEMBLY_INFO.Belem = Belem ;
    ASSEMBLY_INFO.DOFsKEEP = DOFsKEEP ;
    ASSEMBLY_INFO.ndim = ndim ;
    ASSEMBLY_INFO.nstrain = nstrain ;
    ASSEMBLY_INFO.ngaus = ngaus ;
    ASSEMBLY_INFO.wST = wST ;
    ASSEMBLY_INFO.CN = CN ; 
    
else
    ASSEMBLY_INFO = [] ;
end


 
