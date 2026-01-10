function [K,TRANSF_COORD] = ComputeK_EIFEnv(COOR,CN,Bmat,WEIGHTS, PROPMAT,MaterialType,TRANSF_COORD) ;
%%%% STIFFNESS MATRIX, EIFE method, nonvectorized version 
% JAHO, 12-March-2023
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/01_HOMOG/README_HOMOG.mlx
if nargin == 0
    load('tmp.mat')
end
nnode = size(COOR,1); ndim = size(COOR,2); nelem = size(CN,1); nnodeE = size(CN,2) ;
% nstrain = size(celasglo,1) ;
% Shape function routines (for calculating shape functions and derivatives)
% TypeIntegrand = 'K';
% [weig,posgp,shapef,dershapef] = ComputeElementShapeFun(TypeElement,nnodeE,TypeIntegrand) ;
% Assembly of matrix K
% ----------------
% Number of nonzero elements ? We have nelem elements, and each element is
% a (nnodeE*ndim) x (nnodeE*ndim) matrix. Therefore
nzeros = (nnodeE*ndim)^2*nelem ;
K = sparse([],[],[],nnode*ndim,nnode*ndim,nzeros) ;
%Vrot = cell(nelem,1) ;
%dbstop('22')
for e = 1:nelem 
  
    IndexDomainLOC =  TRANSF_COORD{e}.IndexParentDomain ;   
    EIFEoper = PROPMAT(MaterialType(e)).EIFE_prop(IndexDomainLOC) ;   
   
     Ke = ComputeKeEIFEnv(EIFEoper,Bmat{e},WEIGHTS{e}) ;
            
    for anod=1:nnodeE
        a = Nod2DOF(anod,ndim) ;
        for bnod= 1:nnodeE
            b = Nod2DOF(bnod,ndim) ;
            Anod = CN(e,anod) ;  A = Nod2DOF(Anod,ndim) ;
            Bnod = CN(e,bnod) ;  B = Nod2DOF(Bnod,ndim) ;
            %%%%%
            % dbstop('35')
            % Method 1
            K(A,B) = K(A,B) + Ke(a,b) ;
            % Method 2 :  % Accelerate assembly using sparse structure
            % ------------------------------------------
            %             ss = Ke(a,b); ss = ss(:);
            %             ii = repmat(A,size(A,1),1) ;
            %             jj = repmat(B',size(B,1),1) ;  jj = jj(:);
            %    K = K + sparse(ii,jj,ss,nnode*ndim,nnode*ndim,ndim^2) ;
        end
    end
    
    if mod(e,1)==0
        disp(['e=',num2str(e)])
    end
end


