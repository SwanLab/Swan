function GeneralizedForces = generalized_forces_onlyslices(DATA_REFMESH,DATAROM,rDEF,rRB,ndim)

V = DATAROM.BasisINT ; 
 nnode = size(MESH1D.COOR,1) ; % Number of 1D nodes (interfaces)
nelem = size(MESH1D.CN,1)  ; % Number of elements (slices)
nnodeE = 2; % Number of nodes per element (number of interfaces per element)


M = DATA_REFMESH.M ;
f1 = DATAROM.f1 ;
f2 = DATAROM.f2 ;
Tdef_1 = V'*DATAROM.BasisRdef(DATAROM.f1,:) ;
Trb_1 = V'*M(f1,f1)*DATA_REFMESH.BasisUrb(DATAROM.f1,:) ;
Tdef_2 = V'*DATAROM.BasisRdef(DATAROM.f2,:) ;
Trb_2 = V'*M(f2,f2)*DATA_REFMESH.BasisUrb(DATAROM.f2,:) ;

GeneralizedForces = zeros(6,nnode) ;

% This function should be adapted to cover any connectivity topology
rRB = cell2mat(rRB) ;
rRB = reshape(rRB,size(Trb_1,2),[]) ;
rDEF = cell2mat(rDEF) ;
rDEF = reshape(rDEF,size(Tdef_1,2),[]) ;

GeneralizedForces(:,1:nelem) = Tdef_1*rDEF + Trb_1*rRB ;

GeneralizedForces(:,end) = -Tdef_2*rDEF(:,end) - Trb_2*rRB(:,end) ;

