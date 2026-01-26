function GF_FACE= ForcesMoments_PLATE(DATAROM,MESH2D,rDEF,rRB,DATA_REFMESH,ndimINTF) ;


if nargin ==0
    load('tmp0.mat')
end



nnode = size(MESH2D.COOR,1) ; % Number of 1D nodes (interfaces)
nelem = size(MESH2D.CN,1)  ; % Number of elements (slices)
nnodeE = size(MESH2D.CN,2); % Number of nodes per element (number of interfaces per element)

nentities = length(DATAROM) ;
nmodes = 6 ;
GeneralizedForces = cell(1,nnodeE) ;
GeneralizedForces(:) = {zeros(nmodes,nelem)} ;  % Generalized forces
% --------------------------------------------------------
rRB = cell2mat(rRB) ;
rRB = reshape(rRB,nmodes,[]) ;
nmodesR = size(rDEF,1) ;
%rDEF = cell2mat(rDEF) ;
%rDEF = reshape(rDEF,nmodesR,[]) ;

for ientity = 1:nentities
    ELEMS = find(MESH2D.MaterialType == ientity) ;  %
    %%%%% We shall distinguish between forces and moments applied on faces
    %%%%% A,B (parallel to plane x = 0), and those applied on faces C,D
    % Let us start by computing the forces applied to faces A,B
    
    for iface = 1:nnodeE
        GeneralizedForces{iface}(:,ELEMS) = ComputeGeneralizeForce(DATA_REFMESH,ientity,iface,rDEF,rRB,ELEMS,DATAROM,nmodes) ;
    end
end

%%%% Now we are going to define "average" generalized forces 
Lx = MESH2D.xMAX(1)-MESH2D.xMIN(1) ;
Ly = MESH2D.xMAX(2)-MESH2D.xMIN(2) ;

FX =  0.5/Ly*(GeneralizedForces{1} + GeneralizedForces{3}) ;  
GF_FACE.FORCES_PLANEX.VALUE = FX(1:3,:) ; 
GF_FACE.FORCES_PLANEX.LEGEND = {'Nxx','Nxy','Nxz'} ; 
GF_FACE.FORCES_PLANEX.NAME = 'Forces plane X' ; 
GF_FACE.MOMENTS_PLANEX.VALUE = FX(4:6,:) ; 
GF_FACE.MOMENTS_PLANEX.LEGEND = {'Mxx','Mxy','Mxz'} ; 
GF_FACE.MOMENTS_PLANEX.NAME = 'Moments plane X' ;
FY =  0.5/Lx*(GeneralizedForces{2} + GeneralizedForces{4}) ;  
GF_FACE.FORCES_PLANEY.VALUE = FY(1:3,:) ; 
GF_FACE.FORCES_PLANEY.LEGEND = {'Nyx','Nyy','Nyz'} ; 
GF_FACE.FORCES_PLANEY.NAME = 'Forces plane Y' ;

GF_FACE.MOMENTS_PLANEY.VALUE = FY(1:3,:) ; 
GF_FACE.MOMENTS_PLANEY.LEGEND = {'Myx','Myy','Myz'} ; 
GF_FACE.MOMENTS_PLANEY.NAME = 'Moments plane Y' ;


end


function GeneralizedForces_loc = ComputeGeneralizeForce(DATA_REFMESH,ientity,iface,rDEF,rRB,ELEMS,DATAROM,nmodes)

M = DATA_REFMESH{ientity}.M ;
V = DATAROM{ientity}.BasisIntRB{iface} ;
V = V(:,1:nmodes) ;
f = DATAROM{ientity}.fI{iface} ;
Tdef = V'*DATAROM{ientity}.BasisRdef(f,:) ;
Trb = V'*M(f,f)*DATA_REFMESH{ientity}.BasisUrb(f,:) ;
GeneralizedForces_loc = Tdef*rDEF(:,ELEMS) + Trb*rRB(:,ELEMS) ;


end


