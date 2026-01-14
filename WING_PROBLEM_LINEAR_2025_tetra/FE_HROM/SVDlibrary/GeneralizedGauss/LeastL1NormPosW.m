function delta_q =  LeastL1NormPosW(DATALOC,D,Fk,w,m) % 
% See LeastL1NormPosW_aux.mlx
if nargin == 0
    load('tmp1.mat')
end 

mPOINTS = m; 


% SVD on the Jacobian matrix 
 DATALOCSVD.RELATIVE_SVD = 1;  
 TOL = 1e-10 ; 
    [U,S,V] = SVDT(D,TOL,DATALOCSVD) ;     
 U = bsxfun(@times,U',1./S)'  ;
 % Independent term appearing in P*y =Fu
 Fu = U'*Fk ; 
 % Number of variables
 m = size(D,2) ; 
 % Matrix E appearing in the objective function  E^T*y 
 E= [zeros(m,1); ones(m,1)]; 
 E = sparse(E) ; 
 % Matrix P such that P*y = Fu
 r = size(V,2) ; 
 P = [V',zeros(r,m)]; 
 P = sparse(P); 
 % Inequality constraint D*y <= h 
 mPOINTS_nDIM = m-mPOINTS ; 
 C = [zeros(mPOINTS,mPOINTS_nDIM),-eye(mPOINTS)] ; 
  
 
 D=  [eye(m),-eye(m)
      -eye(m),-eye(m)
      C   , zeros(mPOINTS,m)
      zeros(m,m),-eye(m)] ; 
  D = sparse(D) ; 
  
  h = [zeros(m,1); zeros(m,1); w; zeros(m,1)];  
  
  % LINEAR PROGRAMMING 
  
 options = optimoptions('linprog','Algorithm','dual-simplex');
x0 = [];
  lb = []; ub = [];
[y,fval,exitflag,output,lambda] = linprog(E,D,h,P,Fu,lb,ub,x0,options) ; 


if  exitflag ==1
delta_q = y(1:m) ; 
else
     delta_q = V'\Fu ;
end


 
 
        