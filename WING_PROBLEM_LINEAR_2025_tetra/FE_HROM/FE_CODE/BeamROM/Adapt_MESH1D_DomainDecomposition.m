function [  MESH1Ddom] = Adapt_MESH1D_DomainDecomposition(MESH1D,MESH3D)

% Preparing 1D MESH information for the Domain Decomposition problem

% --------------------------
if nargin == 0
    load('tmp1.mat')
end
MESH1Ddom = MESH1D ;

%MESH1Ddom.PROP = MESH1D.PROP ;
%MESH1Ddom.NODES_LINES = MESH1D.NODES_LINES;
%MESH1Ddom.NODES_POINTS =MESH1D.NODES_POINTS ;
%MESH1Ddom.COOR = MESH1D.COOR;
for ientity = 1:length(MESH1Ddom.PROP)
    switch  MESH1Ddom.PROP(ientity).TYPE
        case 'BEAM'
            islice= MESH1Ddom.PROP(ientity).INDEX_ENTITY ;
            MESH1Ddom.PROP(ientity).NameWSmodes = MESH3D.SLICES(islice).NameWSmodes ;
        case 'JOINT'
            ijoint= MESH1Ddom.PROP(ientity).INDEX_ENTITY ;
            MESH1Ddom.PROP(ientity).NameWSmodes = MESH3D.JOINTS(ijoint).NameWSmodes ;
    end
end



nnodeINI = size(MESH1D.COOR,1) ;
% % New connectivity matrix
% ---------------------------
nelemINI = size(MESH1D.CN,1) ;
nnodeE = size(MESH1D.CN,2) ;
nnodeEmax = 10 ; % Maximum number of nodes per joint (normally it is not greater than 4)
CN =  zeros(nelemINI,nnodeEmax);
CN(:,1:nnodeE) = MESH1D.CN  ;
% ------- BEAMS --------------------------------------------
% For beams, we shall switch   the order of the connectivity
% matrix so that  local x axis coincides with the direction
% indicated by the connectivity matrix
% This information was already generated when constructing the
% rotation matrix
%  ---> MESH1D.ORDER_CONNECTIVITIES
% Therefore, we only have to switch the corresponding columns
ORDER_CN = MESH1D.ORDER_CONNECTIVITIES ;
INDEX_CHANGE = find(ORDER_CN == -1) ;
if ~isempty(INDEX_CHANGE)
    CN(INDEX_CHANGE,1) = MESH1D.CN(INDEX_CHANGE,2) ;
    CN(INDEX_CHANGE,2) = MESH1D.CN(INDEX_CHANGE,1) ;
     
end

% NEW VARIABLES
NODES_INTERFACES = ones(nnodeINI,1) ; % Nodes actually coming into play in the physical problem
NNODES = 2*ones(nelemINI,1) ; % Number of nodes per entity (beam = 2)
ISBEAM = ones(nelemINI,1) ; % 1= Beam, 0 = Joint
TYPE_BEAM_JOINT_INDEX = MESH1D.MaterialType ; % Index associated to the element (specified in GID file)
SUBTYPE_BEAM_JOINT_INDEX = zeros(size(MESH1D.MaterialType)) ; % Portions or beams/Joints
SLICE_JOINT_INDEX = MESH1D.TypeOfSlices ; % Index identifying the type of slices or joints  (3D entities)
% Elements to print
% -----------------
Elements2Print = zeros(nelemINI,1) ;
Elements2Print(MESH1D.Elements2Print) =1  ;
% Elements to eliminate (joints)
% ---------------------
Elements2Eliminate = zeros(nelemINI,1) ;

% ------------- JOINTS -------------------------------
% NExt step consists in identifying 1D elements pertaining to JOINTS
% We shall remove these elements except for one of them, in which the end
% nodes of the joint will be stored
% -----------------------------------------------------------------------
for iENTITY  = 1: length(MESH1D.PROP) % Loop over entities (either beam or joints)
    switch MESH1D.PROP(iENTITY).TYPE ; % See whether it is a beam or a joint
        case {'BEAM'}
            ELEMENTS = MESH1D.INFOLINES.ELEMENTS{iENTITY} ;
            for ibeam = 1:length(ELEMENTS)
                ELEMloc  = ELEMENTS{ibeam} ;
                SUBTYPE_BEAM_JOINT_INDEX(ELEMloc) = ibeam ;
            end
            
        case {'JOINT'}
            % Elements pertaining to this type of joints
            ELEMSall = find(MESH1D.MaterialType == iENTITY) ;
            % Nodes pertaining to this type of joints
            NODESall = unique(MESH1D.CN(ELEMSall,:)) ;
            ELEMSkeep = [] ;
            NODESkeep = [] ;
            ISBEAM(ELEMSall) = 0 ; % These elements are JOINTS
            SLICE_JOINT_INDEX(ELEMSall) = MESH1D.PROP(iENTITY).INDEX_ENTITY ; % Index 3D entity
            %%%%%% Information of each joint
            END_ELEMENTS = MESH1D.INFOLINES.END_ELEMENTS{iENTITY} ;
            END_NODES = MESH1D.INFOLINES.END_NODES{iENTITY} ;
            njoints = length(END_NODES) ;
            % Loop over joints
            % -----------------
            for ijoint = 1:njoints
                EndNodesLocal = END_NODES{ijoint} ;
                EndElemLocal = END_ELEMENTS{ijoint} ;
                % We keep one rows of the CN matrix for each joint
                % First element (face 1)
                FIRST_ELEMENT = EndElemLocal(1) ;
                % Therefore
                CN(FIRST_ELEMENT,1:length(EndNodesLocal)) = EndNodesLocal ;
                SUBTYPE_BEAM_JOINT_INDEX(FIRST_ELEMENT) = ijoint ;
                NNODES(FIRST_ELEMENT) = length(EndElemLocal) ;
                % Remaining nodes and elements can be eliminated
                NODESall  = setdiff(NODESall,EndNodesLocal)  ;
                ELEMSall  = setdiff(ELEMSall,FIRST_ELEMENT)  ;
                
            end
            % nodes to remove
            NODES_INTERFACES(NODESall) = 0 ;
            Elements2Eliminate(ELEMSall) = 1;
            
    end
    
    
end


% FINAL CONNECTIVITY MATRIX
% -------------------------


% MAx of nnode
MaxNodeElements  = max(NNODES) ;
% SElected elements
ELEMS_eliminate =  find(Elements2Eliminate==1) ;
ELEMS_SELECTED = setdiff(1:length(Elements2Eliminate),ELEMS_eliminate) ;
MESH1Ddom.ELEMENTS_MESH1D = ELEMS_SELECTED ;
MESH1Ddom.NNODES_elem = NNODES(ELEMS_SELECTED) ;
% ------------------
% NEW CONNECTIVITIES (considering slices as type of super-finite elements )
% ------------------
MESH1Ddom.CN = CN(ELEMS_SELECTED,1:MaxNodeElements) ;
MESH1Ddom.MaterialType = MESH1D.MaterialType(ELEMS_SELECTED) ;
% -------------------------
% List of interface nodes
% -------------------------
MESH1Ddom.NODES_INTERFACES  = unique(find(NODES_INTERFACES==1 )) ;
%
%  First column: 1 == IS beam , 0 is a joint
% Second colums: 3D domain associated to each type (either JOINT of SLICE )
MESH1Ddom.TYPE.ISBEAM = ISBEAM(ELEMS_SELECTED) ;
MESH1Ddom.TYPE.INDEX3D = SLICE_JOINT_INDEX(ELEMS_SELECTED) ;
MESH1Ddom.TYPE.INDEX_BEAM_JOINT = TYPE_BEAM_JOINT_INDEX(ELEMS_SELECTED) ;
MESH1Ddom.TYPE.SUBINDEX_BEAM_JOINT = SUBTYPE_BEAM_JOINT_INDEX(ELEMS_SELECTED) ;


% -----------------------
% Rotation matrix domains and transformation matrices
% -----------------------
ndim = size(MESH1Ddom.COOR,2);
nelem = length(ELEMS_SELECTED) ;
if ~isempty(MESH1Ddom.TRANSFM)
MESH1Ddom.ROTATIONS = zeros(ndim,nelem*ndim) ;
MESH1Ddom.TRANSFM.A = zeros(ndim,nelem*ndim) ;
MESH1Ddom.TRANSFM.D = zeros(ndim,nelem*ndim) ;
MESH1Ddom.TRANSFM.a0 = zeros(ndim,nelem) ;
end
for idim  = 1:ndim
    indexOLD =  ndim*ELEMS_SELECTED-ndim+idim ;
    indexNEW = idim:ndim:nelem*ndim ;
    MESH1Ddom.ROTATIONS(:,indexNEW) =   MESH1D.ROTATIONS(:,indexOLD) ;
    
    if ~isempty(MESH1Ddom.TRANSFM)
    MESH1Ddom.TRANSFM.A(:,indexNEW) =   MESH1D.TRANSFM.A(:,indexOLD) ;
    MESH1Ddom.TRANSFM.D(:,indexNEW) =   MESH1D.TRANSFM.D(:,indexOLD) ;
    
    MESH1Ddom.TRANSFM.a0 =   MESH1D.TRANSFM.a0(:,ELEMS_SELECTED) ;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Coordinate matrix and connectivity matrix to be used in the assembly of
% % the  stiffness matrix
% %
% MESH1dom.CN_plot = MESH1Ddom.CN ;      % Employed for plotting purposes
% MESH1dom.COOR_plot = MESH1D.COOR ;
% % New coordinate and connectivity matrix
% NewNodes = unique(MESH1Ddom.CN(:) ) ;
% MESH1dom.CN = RenumberConnectivities(MESH1Ddom.CN,1:length(NewNodes)) ;  % Renumbering connectivities
% MESH1dom.COOR = MESH1D.COOR(NewNodes,:) ;
%
% MESH1Ddom.NODES_INTERFACES  = 1:length(NewNodes) ;
%
% % Change also POINTS defined for boundary conditions
% MESH1Ddom.NODES_LINES = MESH1D.NODES_LINES;
% MESH1Ddom.NODES_POINTS =MESH1D.NODES_POINTS ;
% for i = 1:length(MESH1D.NODES_LINES)
%
% end
% for i = 1:length(MESH1D.NODES_POINTS)
%
% end









% IT ONLY REMAINS TO DEFINE THE ROTATION MATRICES FOR THE INTERFACE NODES
% ------------------------------------------------------
MESH1Ddom.ROTATIONSint = zeros(ndim,size(MESH1Ddom.COOR,1)*ndim) ;
% By default, this rotation matrices are taken as those of the domains
% whose face F1 is in the node under consideration
% ---- With this criterion, the rotation matrices of certain nodes shared
% by joints will be undefined. In such cases, the rotation matri qx
% associated to face F2 will be taken.

BEAMelements = find(MESH1Ddom.TYPE.ISBEAM == 1) ;
VISITED_NODE = zeros(size(MESH1Ddom.COOR,1),1) ;
for ielemLOC = 1:length(BEAMelements)
    ielem = BEAMelements(ielemLOC) ;
    % Rotation matrix
    ifin = ielem*ndim ; iini = ifin-ndim+1;
    R = MESH1Ddom.ROTATIONS(:,iini:ifin) ;
    % First node of the beam element
    inode = MESH1Ddom.CN(ielem,1) ;
    % Matrix R is assigned to this node (interface node)
    ifin = inode*ndim ; iini = ifin-ndim +1 ;
    MESH1Ddom.ROTATIONSint(:,iini:ifin) = R ;
    VISITED_NODE(inode) = 1;
end

% Check visited nodes
VISITED_NODE = find(VISITED_NODE) ;
% Find out which nodes has not been visited ---typically end nodes
NON_VISITED_NODES = setdiff(MESH1Ddom.NODES_INTERFACES ,VISITED_NODE) ;

%%
for inodeLOC = 1:length(NON_VISITED_NODES)
    inode = NON_VISITED_NODES(inodeLOC) ; %
    [ielemGLO,dummy  ] = find(MESH1Ddom.CN == inode) ;  % Elements sharing this node
    if length(ielemGLO) == 1
        % This node is an end-node (boundary node)
        OTHER_NODE = MESH1Ddom.CN(ielemGLO,:) ;
        OTHER_NODE = setdiff(OTHER_NODE,inode) ;
    end
    % Loop over ielemGLO  to see which element is a beam element
    ielemLOC = 1;
    FOUND = 0 ;
    while   ielemLOC <=length(ielemGLO)
        ielem = ielemGLO(ielemLOC) ;
        if MESH1Ddom.TYPE.ISBEAM(ielem)
            FOUND = 1;
            % First we determine the rotation matrix of the element
            ifin = ielem*ndim ; iini = ifin-ndim+1;
            R = MESH1Ddom.ROTATIONS(:,iini:ifin) ;
            % R is the rotation matrix whose first column is a vector
            % normal to FACE F1. But we are interested in the rotation
            % matrix of FACE F2. The connection between both rotation matrices is
            % given by the rotation matrix derived from the normals between
            % face F1 and face F2. The 3D slice associated to this element has the index:
            INDEX_SLICE  = MESH1Ddom.TYPE.INDEX3D(ielem);
            % Relative rotation
            DeltaR = MESH3D.SLICES(INDEX_SLICE).DATA3D.ROTATIONf1f2 ;
            % This matrix is computed in SliceNodesFacesIdentification. The
            % first vector of this matrix is the normal to face f2
            % expressed in the coordinates attached to face f1.  Therefore,
            % we can write DeltaR = [DeltaR]_{F1,F2}, that is, it is a
            % matrix that takes a vector expressed in the coordinates of
            % face F2 and gives the coordinates of the same vector on the
            % basis of F1.
            % On the other hand, R = [R]_{G,F1}, that is, it takes vectors
            % expressed in F1 and returns the coordinates of the vector on
            % the global system. Here we are interested in [R]_{G,F2},
            % thus, [R]_{G,F2} = [R]_{G,F1}*[DeltaR]_{F1,F2}
            % Thus
            R_f2 = R*DeltaR ;
            ifin = inode*ndim ; iini = ifin-ndim +1 ;
            MESH1Ddom.ROTATIONSint(:,iini:ifin) = R_f2;
            
            
            break
        end
        ielemLOC = ielemLOC + +1;
    end
    if FOUND == 0
        % error(['It was not possible to define the interface rotation matrix of node =',num2str(inode)])
        % This means that this node is an "end" node. In this case, we take
        % as rotation matrix that of the other node of the domain (OTHER_NODE)
        ifin = inode*ndim ; iini = ifin-ndim +1 ;
        ifinOTHER = OTHER_NODE*ndim ; iiniOTHER = ifinOTHER-ndim +1 ;
        
        MESH1Ddom.ROTATIONSint(:,iini:ifin) = MESH1Ddom.ROTATIONSint(:,iiniOTHER:ifinOTHER) ;
    end
    
end


for ientity = 1:length(MESH1D.PROP)
    MESH1Ddom.PROP(ientity).xLOCAL = 1 ;  % This is already included in the connectivity matrix 
end


end

