function P = AssemblyPinterfacesRVE_PLATE(DATAROM,MESH2D,DATAIN,FORCES,DATA_REFMESH,ndim)

if nargin == 0
    load('tmp0.mat')
end

nnode = size(MESH2D.COOR,1) ; % Number of 1D nodes (interfaces)
nelem = size(MESH2D.CN,1)  ; % Number of elements (slices)
nnodeE = length(DATAROM{1}.BasisIntRB); % Number of nodes per element (number of interfaces per element)


%%% Matrix that maps forces and moments defined over faces to forces and
%%% moments over corners
[Js1,Js2] = MapsFacesCorners(DATA_REFMESH) ; 




% STEP 1
P = zeros(nnode*ndim,1) ;

%
NODES_LINESall = MESH2D.NODES_LINESall ;  % Both midside and corner points
nlines = min(length(NODES_LINESall),length(FORCES.LINE)) ;

for iline = 1:nlines  % loop  over the set of lines defined by the user (through GID)
  %  FORCESloc = FORCES.LINE{iline} ;   % Generalized forces, 6-by-1, applied to this line
   [xC,yC,L_total,CNb] = CoordLineElement(NODES_LINESall,MESH2D,iline) ; 
    
    for ielem = 1:size(CNb,1)        
        % Choosing transformation matrix, and computing the associated DOFS
        % to the given line 
        [J,DOFS,LENGTH] = BoundaryCond_DOFS_and_J_M(MESH2D,CNb,ielem,Js1,Js2,ndim,xC,yC) ;
        % Computing FORCESloc (forces with respect to the centroid of the boundary element)
        FORCESloc = FORCES.LINE{iline} ; 
        %% ASSEMBLY        
        Pe = J*(FORCESloc*LENGTH/L_total) ;        
        P(DOFS) = P(DOFS) + Pe ;   
    end
    
    
    
    
    %     DOFS = small2large(NODES_LINES{iline},ndim) ;
    %
    %     nnodes = length(NODES_LINES{iline}) ;
    %     FORCES_rep = repmat(FORCESloc,nnodes,1) ;
    %     P(DOFS) = FORCES_rep ;
end




end




