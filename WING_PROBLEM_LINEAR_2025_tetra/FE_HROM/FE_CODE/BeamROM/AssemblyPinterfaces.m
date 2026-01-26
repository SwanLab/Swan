function P = AssemblyPinterfaces(DATAROM,MESH1D,DATAIN,FORCES,DATA_REFMESH,ndim)

if nargin == 0
    load('tmp4.mat')
end



nnode = size(MESH1D.COOR,1) ; % Number of 1D nodes (interfaces)
nelem = size(MESH1D.CN,1)  ; % Number of elements (slices)
nnodeE = 2; % Number of nodes per element (number of interfaces per element)

% Node left/right end
% -------------
MESH1D= DefaultField(MESH1D,'LEFT_END_NODE',[]) ; 
if  ~isempty(MESH1D.LEFT_END_NODE )
    NODE{1} = MESH1D.LEFT_END_NODE ;
    NODE{2} = MESH1D.RIGHT_END_NODE ;
    FORCES_loc{1} = FORCES.LEFT_END;
    FORCES_loc{2} = FORCES.RIGHT_END;
    
else
    NODE{1} = MESH1D.NODES_POINTS{1} ;
    NODE{2} = MESH1D.NODES_POINTS{2} ;
    FORCES_loc{1} = FORCES.NODES{1};
    FORCES_loc{2} = FORCES.NODES{2};
end
P = zeros(nnode*ndim,1) ;


for ifaces = 1:length(NODE)
    DOFS{ifaces} = small2large(NODE{ifaces},ndim) ;  
    nforces =length(FORCES_loc{ifaces}) ;
    P(DOFS{ifaces}(1:nforces)) = FORCES_loc{ifaces} ;    
end