function [GeneralizedForces,RotationMatrixFace2 ]= generalized_forces_jslices(DATA_REFMESH,DATAROM,rDEF,rRB,MESH1D,DATAIN)

if nargin == 0
    load('tmp2.mat')
end


nnode = size(MESH1D.COOR,1) ; % Number of 1D nodes (interfaces)
nelem = size(MESH1D.CN,1)  ; % Number of elements (slices)
nnodeE = 2; % Number of nodes per element (number of interfaces per element)

nentities = length(DATAROM) ;
nmodesRB = size(rRB{1},1) ;

DATAROM{1} = DefaultField(DATAROM{1},'WarpingModesInterface',[]) ; 

if ~isempty(DATAROM{1}.WarpingModesInterface)
    nmodes = nmodesRB+1; % Include warping
else
    nmodes = nmodesRB ;
end

GeneralizedForces = zeros(nmodes,nnode) ;  % Generalized forces

% --------------------------------------------------------
% Only valid for straight, mixed beams (slices and joints)
rRB = cell2mat(rRB) ;
rRB = reshape(rRB,nmodesRB,[]) ;
nmodesR = size(rDEF,1) ;
%rDEF = cell2mat(rDEF) ;
%rDEF = reshape(rDEF,nmodesR,[]) ;

for ientity = 1:nentities
    
   
    ELEMS = find(MESH1D.MaterialType == ientity) ;  % Elements featuring this type of slice/JOINT
  %  NODES_FACE1 = MESH1D.CN(ELEMS,iface) ; 
    [Tdef,Trb,RotationMatrixFace2] = ComputeOperatorsGeneralizedForces(MESH1D,DATAROM,ientity,nmodes,DATA_REFMESH,ELEMS) ; 
    
     iface = 2;
    for ielem = 1:length(ELEMS)
        nodeloc  =  MESH1D.CN(ELEMS(ielem),iface) ;
        Force_1 =  Tdef{iface}*rDEF(:,ELEMS(ielem)) - Trb{iface}*rRB(:,ELEMS(ielem)) ;
        %  Force_2 =  Tdef_2*rDEF(:,ELEMS(ielem)) + Trb_2*rRB(:,ELEMS(ielem)) ;  %% MINus signus is necessary
        GeneralizedForces(:,nodeloc) =Force_1 ;
    end
    
end






% Adaption of old code 
MESH1D  =DefaultField(MESH1D,'RIGHT_END_NODE',[]) ; 
if ~isempty(MESH1D.RIGHT_END_NODE)
    MESH1D.NODES_POINTS{1} = MESH1D.RIGHT_END_NODE ; ;
     MESH1D.NODES_POINTS{2} = MESH1D.LEFT_END_NODE ; ;
end

%   
ENDNODE_REF =  MESH1D.NODES_POINTS{1}  ; 
[ELEMS jjj]= find(MESH1D.CN == ENDNODE_REF) ;
if jjj==1  
    % This means that it is the other node
    ENDNODE  =   MESH1D.NODES_POINTS{1} ;
else
    ENDNODE  =   MESH1D.NODES_POINTS{2} ;
end

%% FACE 1
[ELEMS jjj]= find(MESH1D.CN == ENDNODE) ;
ientity = MESH1D.MaterialType(ELEMS) ;
iface=1 ; 

[Tdef,Trb] = ComputeOperatorsGeneralizedForces(MESH1D,DATAROM,ientity,nmodes,DATA_REFMESH,ELEMS) ; 
GeneralizedForces(:,ENDNODE) = -Tdef{iface}*rDEF(:,ELEMS) + Trb{iface}*rRB(:,ELEMS) ;


