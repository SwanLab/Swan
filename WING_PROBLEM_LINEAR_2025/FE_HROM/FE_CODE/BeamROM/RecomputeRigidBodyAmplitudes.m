function  qRB =  RecomputeRigidBodyAmplitudes(qDEF,DATA_REFMESH,DATAROM,a) 

if nargin == 0 
    load('tmp.mat')
end
DATAROM = DATAROM{1} ; 
DATA_REFMESH = DATA_REFMESH{1} ; 
%%% RIGID BODY AMPLITUDES 
f1 = DATAROM.f1 ;  % DOFs face 1
f2 = DATAROM.f2 ;  % DOFs face 2
f = [f1;f2] ; 
R  = DATA_REFMESH.BasisUrb  ; % Rigid body modes faces 1 and 2 
M = DATA_REFMESH.GeometricMassMatrixInterface ; 
M = M{1} ; %ç
ndim = 3; 
V= DATAROM.BasisInt ; 
% M3D = zeros(ndim*size(M)) ; 
% for idim = 1:ndim 
%     M3D(idim:ndim:end,idim:ndim:end) = M ; 
% end

% Rf1*M*Rf1 
nmodesRB = size(R,2) ; 
Rf1 = R(f1,:) ; 
Rf2 = R(f2,:) ; 
Rf1_f1 = zeros(nmodesRB,nmodesRB) ; 
Rf1_f2 = zeros(nmodesRB,nmodesRB) ; 
Rf2_f1 = zeros(nmodesRB,nmodesRB) ;
Rf2_f2 = zeros(nmodesRB,nmodesRB) ;

Rf1_m = zeros(size(Rf1)) ; 
Rf2_m = zeros(size(Rf2)) ; 


for idim = 1:ndim
    Rf1_f1 = Rf1_f1+ Rf1(idim:ndim:end,:)'*M*Rf1(idim:ndim:end,:) ; 
    Rf1_f2 = Rf1_f2+ Rf1(idim:ndim:end,:)'*M*Rf2(idim:ndim:end,:) ; 
    Rf2_f1 =Rf2_f1 +  Rf2(idim:ndim:end,:)'*M*Rf1(idim:ndim:end,:) ; 
    Rf2_f2 = Rf2_f2 + Rf2(idim:ndim:end,:)'*M*Rf2(idim:ndim:end,:) ; 
    Rf1_m(idim:ndim:end,:) =   M*Rf1(idim:ndim:end,:) ;
    Rf2_m(idim:ndim:end,:) =   M*Rf2(idim:ndim:end,:) ;
end

% Computing deformational displacements face f1
%qDEF = cell2mat(qDEF') ; 
dDEF_f1 = DATAROM.BasisUdef(f1,:)*qDEF ; 
dDEF_f2 = DATAROM.BasisUdef(f2,:)*qDEF ; 

DeltaDEF = dDEF_f2(:,1:end-1) - dDEF_f1(:,2:end) ;
nDOM = size(qDEF,2) ; 
b1 = zeros(nmodesRB,nDOM) ; 
b2 = b1 ; 

% Vector b
% ----------
b1(:,2:end) = -Rf1_m'*DeltaDEF ; 
b2(:,1:end-1) = Rf2_m'*DeltaDEF ; 
b = b1+b2 ; 
b = b(:) ; 

nmodesINT = size(V,2) ; 
a1 = a(1:nmodesINT) ; 
b(1:nmodesRB) = b(1:nmodesRB) - Rf1_m'*(V*a1 -dDEF_f1(:,1) ) ; 
aEND = a(end-nmodesINT+1:end) ; 
b(end-nmodesRB+1:end) = b(end-nmodesRB+1:end) - Rf2_m'*(V*aEND -dDEF_f2(:,end) ) ; 


% Matrix P 
% --------
 
% Main diagonal 
% ---------------
R_f = Rf1_f1 +Rf2_f2 ; 
DIAG1 =cell(nDOM,1) ; 
DIAG1(:) ={sparse(R_f)} ; 
DIAG1 = blkdiag(DIAG1{:}) ;

% Upper diagonal 
% ----------------
DIAGup =cell(nDOM-1,1) ; 
DIAGup(:) ={-sparse(Rf2_f1)} ; 
DIAGup = blkdiag(DIAGup{:}) ;

% Lower diagonal 
% ----------------
DIAGlow =cell(nDOM-1,1) ; 
DIAGlow(:) ={-sparse(Rf1_f2)} ; 
DIAGlow = blkdiag(DIAGlow{:}) ;
 
% P = sparse(length(b),length(b)) ; 
 
 P = DIAG1 ; 
 
 % Upper diagonal
 iniROW = 1; finROW = size(DIAGup,1) ; 
 inicol = nmodesRB+1; fincol = nmodesRB+size(DIAGup,2) ; 
 P(iniROW:finROW,inicol:fincol) =  P(iniROW:finROW,inicol:fincol)  +  DIAGup ; 
  % Lower diagonal
 iniROW = nmodesRB+1;; finROW =  nmodesRB+size(DIAGup,2) ;  
 inicol = 1; fincol = size(DIAGlow,2) ; 
 P(iniROW:finROW,inicol:fincol) =  P(iniROW:finROW,inicol:fincol)  +  DIAGlow ; 
 
 
 % Solving for qRB 
 qRB = -P\b ; 
 
 qRB = reshape(qRB,nmodesRB,[]) ; 
 

 
