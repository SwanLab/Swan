function [rDEF ]= AmplitudeReactions_jbeam(DATAROM,MESH1D,a,fextBEAMr,ndim)

if nargin ==0
    load('tmp4.mat')
end

%V = DATAROM.BasisInt ;  % Interface modes matrix
%ndim = size(V,2) ; % Number of entries for each node
nnode = size(MESH1D.COOR,1) ; % Number of 1D nodes (interfaces)
nelem = size(MESH1D.CN,1)  ; % Number of elements (slices)
nnodeE = 2; % Number of nodes per element (number of interfaces per element)


nentities = length(DATAROM) ; % Number of joints/slices
KT = cell(nnodeE,nentities) ;



for ientities = 1:nentities
    % OLD IMPLEMENTATION (STRAIGHT DOMAINS)
    %     if iscell(DATAROM{ientities}.BasisInt)
    %         V1 = DATAROM{ientities}.BasisInt{1} ;
    %         V2 = DATAROM{ientities}.BasisInt{2} ;
    %     else
    %         V1 = DATAROM{ientities}.BasisInt ;
    %         V2 = V1;
    %     end
    %     T_1 = V1'*DATAROM{ientities}.BasisRdef(DATAROM{ientities}.f1,:) ;
    %     T_2 = V2'*DATAROM{ientities}.BasisRdef(DATAROM{ientities}.f2,:) ;
    %     KT{1,ientities} = DATAROM{ientities}.Kbeam*T_1' ;
    %     KT{2,ientities} = DATAROM{ientities}.Kbeam*T_2' ;
    
    for iface = 1:nnodeE
        KT{iface,ientities} = DATAROM{ientities}.Kbeam*(DATAROM{ientities}.Tcomp{iface})' ;
    end
    
end


 
% ONLY VALID FOR EQUAL NUMBER OF MODES  
% --------------------------------------
nmodes = size(fextBEAMr{1},1) ; 
fextBEAMr = cell2mat(fextBEAMr) ; 
rDEF = -fextBEAMr ;
rDEF  =reshape(rDEF,nmodes,[]) ;
for  ientity = 1:length(DATAROM)
    ELEMS = find(MESH1D.MaterialType == ientity) ;
    CNlocNOD = MESH1D.CN(ELEMS,:) ;
    for inode = 1:nnodeE
        NODE = CNlocNOD(:,inode) ;
        DOFS = Nod2DOF(NODE,ndim) ;
        aDOFS = a(DOFS) ;
        aDOFS = reshape(aDOFS,ndim,[]) ;
        rDEF(:,ELEMS) =  rDEF(:,ELEMS) + KT{inode,ientity}*aDOFS  ;
    end

end

% VARIABLE NUMBER OF MODES
% rDEF = -fextBEAMr ;
% rDEF  =reshape(rDEF,nmodes,[]) ;
% rDEF = cell(size(fextBEAMr));
% for  ientity = 1:length(DATAROM)
%     ELEMS = find(MESH1D.MaterialType == ientity) ;
%     CNlocNOD = MESH1D.CN(ELEMS,:) ;
%     %rDEF(ELEMS) =-fextBEAMr(ELEMS) ;
%     factor = -1 ; 
%     rDEF = IniCellrDEF(rDEF,ELEMS,fextBEAMr(ELEMS),factor)
%     for inode = 1:nnodeE
%         NODE = CNlocNOD(:,inode) ;
%         DOFS = Nod2DOF(NODE,ndim) ;
%         aDOFS = a(DOFS) ;
%         aDOFS = reshape(aDOFS,ndim,[]) ;
%         %%%
%         rADD = KT{inode,ientity}*aDOFS  ; 
%         %
%         % rDEF(:,ELEMS) =  rDEF(:,ELEMS) + rADD ; 
%         rDEF = AddKT(rDEF,ELEMS,rADD) ; 
%     end
%     
% end


end
% 
% function   rDEF = IniCellrDEF(rDEF,ELEMS,VALUE,factor)
% rLOC = factor*cell2mat(VALUE(ELEMS)')   ;
% [nrows ncols]= size(rLOC) ;
% rLOC = mat2cell(rLOC,nrows,ones(1,ncols));
% rDEF(ELEMS) = rLOC' ;
% end
% 
% function   rDEF = AddKT(rDEF,ELEMS,rADD)
% 
% 
% rDEF_loc = cell2mat(rDEF(ELEMS)') ; 
% nmodes  = size(rDEF_loc,1) ;
% rDEF_loc = rDEF_loc + reshape(rADD,nmodes,[]) ; 
% [nrows ncols]= size(rDEF_loc) ;
% rDEF_loc = mat2cell(rDEF_loc,nrows,ones(1,ncols));
% rDEF(ELEMS) = rDEF_loc' ;
% end


