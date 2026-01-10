function [K,Belem,Cglo] = AssemblyMethodStandard(wSTs,Belem,nstrain,nelem,nnodeE,...
    ndim,CN,nnode,celasglo,wST)
%%  Assembly  method in which the stiffness matrix is computed by first determining
% the product B'(e,xi)*C(e,xi)*B(e,xi), and then performing the assembly in
% the standard, vectorized fashion
% Joaquín A. Hernández (jhortega@cimne.upc.edu), 28-Oct-2015
% ---------------------------------------------------------------------------
%dbstop('9')
if nargin == 0
  load('tmp.mat')
end

ngaus = length(wSTs)/nelem ;
disp('Global elasticity matrix  ...')
% Matrix with all elasticity tensors (at all gauss points)
Cglo = DefineElastMatGLO(celasglo,ngaus)  ; % Cglo =[C(e=1,g=1);C(e=1,g=2) ; ...
%%% 
disp('Point-wise product  C*B')
CB = ProducMatrBlock(Cglo,Belem) ; 
disp('Matrix W*(B)^T') 
% Product of the vector of weights times Belem
BelemT = bsxfun(@times,Belem,wST) ; 
% Transponse
BelemT =   TransposeVectorizeGen(BelemT,ngaus*nelem)  ; 
disp('Point-wise product ...  Bt*C*B')
BCB = ProducMatrBlockGen(BelemT,CB,ngaus*nelem) ;
% 
disp('Matrix of element stiffness matrices')
Kelem = ElementKMatrices(nelem,nnodeE,ndim,BCB,ngaus) ; 
disp('Assembly of K ...')
K = AssemblyKGlobal(Kelem,nelem,nnodeE,ndim,CN,nnode) ;
 
 