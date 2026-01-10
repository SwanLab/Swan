function [ndim,DATAIN]= InformationReducedOrderModel(DATAROM,MESH2D,DATAIN)

if nargin == 0
    load('tmp.mat')
end

 

ndim = [] ;
DATAIN.NODES_ENTITIES  = MESH2D.NODES_LINES;  % 25-Jan-2020



for ientity = 1:length(DATAROM)
    
 
    V = DATAROM{ientity}.BasisInt;
    if ~isstruct(V)
        
        %   if iscell(V)
        for iface = 1:length(V)
            ndim(end+1) = size(V{iface},2) ;
        end
    else
        ndim  = V.nDOFsFACE ; 
    end
    %  else
    %    ndim(end+1) = size(V,2) ;
    % end
    nmodesU = size(DATAROM{ientity}.BasisUdef,2) ;
    nmodesR = size(DATAROM{ientity}.BasisRdef,2) ;
    
   %  nmodesS = size(DATAROM{ientity}.BasisS,2) ;
    nDOF = size(DATAROM{ientity}.BasisUdef,1) ;
    npoints = length(DATAROM{ientity}.HROMVAR.setPoints) ;
    disp('*****************************************************+')
    disp(['Unit Cell =',num2str(ientity)])
    disp(['Number of interfaces = ',num2str(length(ndim))]) ;
    disp(['Number of interface modes (DOFs) = ',num2str(ndim(:)')]) ;
    disp(['Number of disp. modes  = ',num2str(nmodesU)]) ;
    disp(['Number of reaction. modes  = ',num2str(nmodesR)]) ;
%       disp(['Number of stress  modes  = ',num2str(nmodesS)]) ;
    disp(['Number of integration points  = ',num2str(npoints)]) ;
    disp('*****************************************************+')
    
    
end

DATAIN.npointsECM = npoints ; 
DATAIN.ndimINTF = ndim; 

% 25-Jan-2020
% ------------------------------------------------------------------------
% Let us construct a table of "Degrees of freedom". For instance, if there
% is just one element with 7-3-7-3, then 
% NODE 1 : 1,2,...7
% NODE 2: 8 9 10
% Node 3: 11,12,....17
% Node 4: 18,19,20


nnode = size(MESH2D.COOR,1) ;
[ndimMAX faceMAX]= max(ndim) ; 
% Approach based on "ghost" DOFs
% Let us define a matrix ndimMAX x nnode  
DOFsNODEmatrix = ones(ndimMAX,nnode) ; 
% Now let us set to zero those entries corresponding to the "ghost" (non-existing) DOFS
ndimX = ndim(1) ; 
ndimY = ndim(2) ; 
% Which are the faces with less modes ?  --> faceMIN 
% -------------------------------------   -----------
if faceMAX == 1
    faceMIN = [2 4] ;
else
    faceMIN = [1 3] ; 
end
% -------------------------------------------------
% DOFs to eliminate   (reprogram this, not efficient)
% -----------------------------------------
NODESmin = unique(MESH2D.CN(:,faceMIN)) ; 
DOFSmin = small2large(NODESmin,ndimMAX) ; 
DOFSmin = reshape(DOFSmin,ndimMAX,[]) ; 
idofADD = ndim(faceMIN(1))+1:ndimMAX ; 
DOFSeliminate = DOFSmin(idofADD,:) ; 

DOFsNODEmatrix(DOFSeliminate) = 0 ; 


TableDOFSnode = cell(1,nnode) ; 
for innode = 1:nnode
    ncomp = sum(DOFsNODEmatrix(:,innode)); 
    if innode ==1 
     TableDOFSnode{innode} =(1:ncomp)' ;
    else
        nACUM = TableDOFSnode{innode-1}(end) ; 
          TableDOFSnode{innode} = nACUM + (1:ncomp)' ; 
    end
end

DATAIN.TableDOFSnode = TableDOFSnode ; 


