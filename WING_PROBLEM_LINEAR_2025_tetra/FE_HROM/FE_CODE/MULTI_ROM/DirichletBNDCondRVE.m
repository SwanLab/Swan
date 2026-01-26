function [DOFr,dR,DOFm,G,NODES_LINES_ROTATIONS] = DirichletBNDCondRVE(DATAROM,MESH2D,DISP,DATAIN,ndim,DOFsKEEP,DATA_REFMESH)

if nargin == 0
    load('tmp2.mat')
end

ndimSP = size(DATAROM{1})

nnode = size(MESH2D.COOR,1) ; % Number of 2D nodes (interfaces)
nelem = size(MESH2D.CN,1)  ; % Number of elements
nnodeE = length(DATAROM{1}.BasisIntRB); % Number of nodes per element (number of interfaces per element)
[~,nRB] = size(DATAROM{1}.BasisIntRB{1});

% ----------
NODES_LINES = MESH2D.NODES_LINES ;  % Nodes assigned to each "LINE" entity
dR = [] ;
DOFr = [] ;
DOFm = [] ;
nlines = min(length(NODES_LINES),length(DISP.LINE)) ;

ndimMAX = max(ndim) ;
Gcell  = {} ;  G = [] ;  
NODES_LINES_ROTATIONS = cell(1,nlines) ; 
for ilines = 1:nlines
    
    if ~isstruct(DATAROM{1}.BasisInt)
        [DOFrLOC,dRloc,Gloc,DOFmLOC,NODES_LINES_ROTATIONS{ilines}] = BCsSingleLine(NODES_LINES,ilines,DISP,ndimMAX,DATAROM{1}.BasisInt,DATAIN,MESH2D,DATA_REFMESH) ;
    else
        % Method, Apr-2019  (ABANDONED !!! )
        [DOFrLOC,dRloc] = BCsSingleLine_kinemCONS(NODES_LINES,ilines,DISP,ndimMAX,DATAROM{1}.BasisInt,MESH2D,nRB) ;
       
    end
    
    
    DOFr = [DOFr; DOFrLOC] ;
    dR = [dR; dRloc] ;
    DOFm = [DOFm; DOFmLOC]; 
    Gcell{ilines} = Gloc ;
end

G = blkdiag(Gcell{:}) ; 

%%%%%%%%%%%%%%%%%
if ~isempty(DOFsKEEP)
    % Total number of DOFS  ---slave  DOFs
    [DOFrSELECT III JJJ]= intersect(DOFsKEEP,DOFr) ;
    DOFr = III ;
    dR = dR(JJJ) ;
    
    %  
    [DOFmSELECT III JJJmast]= intersect(DOFsKEEP,DOFm) ;
    DOFm = III ;
    
    if ~isempty(G)
    
    G = G(JJJ,JJJmast) ; 
    
    end
    
end

%%%%%%%%%%%%%%%%%%
% 28-JAN-2020, YOU CAN CHECK IF THE DOFS ARE CORRECTLY PRESCRIBED BY
% CONSULTING THE "tABLE OF dofS" -- > For instance
% cell2mat(DATAIN.TableDOFSnode(NODES_LINES{4}))



[DOFr,III] = sort(DOFr) ;
dR = dR(III) ;

[DOFm,JJJ] = sort(DOFm) ;

if ~isempty(G)
G = G(III,JJJ) ; 
end

% DOFl = 1:ndim*nnode ;
% DOFl(DOFr) = [] ;




end
