function Fdis_i = FdisCOMP(COOR,CNb,TypeElementB,tracB,idim) ;
nnode = size(COOR,1);  nelemB = size(CNb,1); nnodeEb = size(CNb,2) ;     ndim = size(COOR,2);
TypeIntegrand = 'RHS';
[weig,posgp,shapef,dershapef] = ComputeElementShapeFun(TypeElementB,nnodeEb,TypeIntegrand) ; 
Fdis_i = zeros(nnode*ndim,1) ; 
for e = 1:nelemB   
    CNloc = CNb(e,:) ;     % Nodes of element "e"   
    tracBe = tracB(e,:)' ; tracBe = tracBe(:) ;   % Nodal values of the distributed load at element "e"    
    Xe = COOR(CNloc,:)' ;% Coordinates of the nodes of element "e"
    Fdis_i_e = FdisElem(tracBe,weig,shapef,dershapef,Xe) ;   % Computation of element traction vector  
    for a=1:nnodeEb 
            Anod = CNb(e,a) ; 
            A = (Anod-1)*ndim +idim ; 
            Fdis_i(A) = Fdis_i(A) + Fdis_i_e(a) ;    
    end
end