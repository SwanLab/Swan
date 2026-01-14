function Fb = ComputeFb(COOR,CN,TypeElement, fNOD) ;
% This subroutine   returns the  body force    contribution (Fb)
%to the % global external force vector. Inputs:   COOR, CN,
%TypeElement, fNOD (nnode*ndim x 1):   Body force function at the nodes  of the mesh.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%dbstop('5')
if nargin==0
    load('tmp2.mat')
end
nnode = size(COOR,1); ndim = size(COOR,2); nelem = size(CN,1); nnodeE = size(CN,2) ;
TypeIntegrand = 'RHS' ;
[weig,posgp,shapef,dershapef] = ComputeElementShapeFun(TypeElement,nnodeE,TypeIntegrand) ;
Fb = zeros(nnode*ndim,1) ;
for e = 1:nelem
    CNlocNOD = CN(e,:) ;
    CNloc = Nod2DOF(CNlocNOD,ndim) ;
    % Body force function evaluated at the nodes of element "e"
    fe = fNOD(CNloc) ;
    % Coordinates of the nodes of element "e"
    Xe = COOR(CNlocNOD,:)' ;
    % Computation of elemental body force  vector
    Fbe = ComputeFbeVector(fe,weig,shapef,dershapef,Xe) ;
    for anod=1:nnodeE
        a = Nod2DOF(anod,ndim) ;
        Anod = CN(e,anod) ; A = Nod2DOF(Anod,ndim) ;
        Fb(A) = Fb(A) + Fbe(a) ;
    end
end