function K = ComputeK(COOR,CN,TypeElement, celasglo) ;
%%%%
% This subroutine   returns the global stiffness matrix K (ndim*nnode x ndim*nnode)
% Inputs:   COOR: Coordinate matrix (nnode x ndim), % CN: Connectivity matrix (nelem x nnodeE),
% TypeElement: Type of finite element (quadrilateral,...),  celasglo (nstrain x nstrain x nelem) 
% Array of elasticity matrices
% Dimensions of the problem
if nargin == 0
    load('tmp1.mat')
end
nnode = size(COOR,1); ndim = size(COOR,2); nelem = size(CN,1); nnodeE = size(CN,2) ;  
% nstrain = size(celasglo,1) ;
% Shape function routines (for calculating shape functions and derivatives)
TypeIntegrand = 'K';
[weig,posgp,shapef,dershapef] = ComputeElementShapeFun(TypeElement,nnodeE,TypeIntegrand) ;
% Assembly of matrix K
% ----------------
K = sparse(nnode*ndim,nnode*ndim) ;
for e = 1:nelem
    celas = celasglo(:,:,e) ;  % Stiffness matrix of element "e"
    CNloc = CN(e,:) ;   % Coordinates of the nodes of element "e"
    Xe = COOR(CNloc,:)' ;     % Computation of elemental stiffness matrix
    Ke = ComputeKeMatrix(celas,weig,dershapef,Xe) ;
    for anod=1:nnodeE
        a = Nod2DOF(anod,ndim) ;
        for bnod= 1:nnodeE
            b = Nod2DOF(bnod,ndim) ;
            Anod = CN(e,anod) ;  A = Nod2DOF(Anod,ndim) ;
            Bnod = CN(e,bnod) ;  B = Nod2DOF(Bnod,ndim) ;
            %%%%%
            K(A,B) = K(A,B) + Ke(a,b) ;
        end
    end
end


