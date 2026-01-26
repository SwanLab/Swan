function P = AssemblyPinterfacesRVE(DATAROM,MESH2D,DATAIN,FORCES,DATA_REFMESH,ndimALL,DOFsKEEP)

if nargin == 0
    load('tmp2.mat')
end


nnode = size(MESH2D.COOR,1) ; % Number of 1D nodes (interfaces)
nelem = size(MESH2D.CN,1)  ; % Number of elements (slices)
nnodeE = length(DATAROM{1}.BasisIntRB); % Number of nodes per element (number of interfaces per element)

% ----------------------
ndim = max(ndimALL) ;



NODES_LINES = MESH2D.NODES_LINES ;
nlines = min(length(NODES_LINES),length(FORCES.LINE)) ;


if isstruct(DATAROM{1}.BasisInt)
    % New method, 27-Apr-2019, Kinem. constraint. 
    % -----------------------
    
    P = BoundaryIntefaceForcesLocal_kinconst(nlines,NODES_LINES,ndim,nnode,DATAROM,MESH2D,DOFsKEEP,FORCES);
   
    
    
else
    
    P = zeros(nnode*ndim,1) ;

    for iline = 1:nlines
        DOFS = small2large(NODES_LINES{iline},ndim) ;
        FORCESloc = zeros(ndim,1) ;
        
         FORCESloc(1:length(FORCES.LINE{iline}))  =     FORCES.LINE{iline} ;
        nnodes = length(NODES_LINES{iline}) ;
        FORCES_rep = repmat(FORCESloc,nnodes,1) ;
        P(DOFS) = FORCES_rep ;
    end
    
    
    if ~isempty(DOFsKEEP)
        P = P(DOFsKEEP) ;
    end
    
    
end

