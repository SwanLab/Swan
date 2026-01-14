function M = ComputeM(COOR,CN,TypeElement, densglo) ;
%%%%
% This subroutine   returns the global mass matrix M (ndim*nnode x ndim*nnode)
% Inputs:   COOR: Coordinate matrix (nnode x ndim), % CN: Connectivity matrix (nelem x nnodeE), 
% TypeElement: Type of finite element (quadrilateral,...),  celasglo (nstrain x nstrain x nelem)  
% densglo: Densities
%dbstop('8')
if nargin == 0
    load('tmp1.mat')
end
nnode = size(COOR,1); ndim = size(COOR,2); nelem = size(CN,1); nnodeE = size(CN,2) ;  
% nstrain = size(celasglo,1) ;
% Shape function routines (for calculating shape functions and derivatives)
TypeIntegrand = 'RHS';
[weig,posgp,shapef,dershapef] = ComputeElementShapeFun(TypeElement,nnodeE,TypeIntegrand) ;
% Assembly of matrix K
% ---------------- 
% Number of nonzero elements ? We have nelem elements, and each element is
% a (nnodeE*ndim) x (nnodeE*ndim) matrix. Therefore
nzeros = (nnodeE*ndim)^2*nelem ; 
M = sparse([],[],[],nnode*ndim,nnode*ndim,nzeros) ;
for e = 1:nelem
    dens = densglo(e) ;  % Stiffness matrix of element "e"
    CNloc = CN(e,:) ;   % Coordinates of the nodes of element "e"
    Xe = COOR(CNloc,:)' ;     % Computation of elemental stiffness matrix
    Me = ComputeMeMatrix(dens,weig,dershapef,Xe,shapef) ;
    for anod=1:nnodeE
        a = Nod2DOF(anod,ndim) ;
        for bnod= 1:nnodeE
            b = Nod2DOF(bnod,ndim) ;
            Anod = CN(e,anod) ;  A = Nod2DOF(Anod,ndim) ;
            Bnod = CN(e,bnod) ;  B = Nod2DOF(Bnod,ndim) ;
            %%%%%         
            M(A,B) = M(A,B) + Me(a,b) ;            
        end
    end
    
    if mod(e,10)==0
        disp(['e=',num2str(e)])
    end
end


