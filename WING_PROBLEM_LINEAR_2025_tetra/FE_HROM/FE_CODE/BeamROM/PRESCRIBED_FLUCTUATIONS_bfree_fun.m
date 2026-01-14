function [G,uBAR,DOFr,DOFm,AREA,R] = PRESCRIBED_FLUCTUATIONS_bfree_fun(a_A,a_B,DOMAINVAR,COOR,CONNECTb,TypeElementB,...
    DATA)

if nargin == 0
    load('tmp1.mat')
end

load(DATA.NameWS_bending_displacements,'dBENDING')  ;




ndim = 3;
%%%% FACE 1
% ----------
iface=1 ;
[DOFA,AREA,R,M,sA,pA,UbA,IbA]...
    =  FLuctBC_onefaceBFREE(iface,DOMAINVAR,COOR,CONNECTb,TypeElementB,dBENDING,...
    a_A,DATA) ;
a_A = cell2mat(a_A) ;
%%%% FACE 2
% ----------
iface=2 ;
nodesfB = DOMAINVAR.NODES_faces12{end,iface} ;
DOFB = small2large(nodesfB,ndim) ;   % DOFS face 2

ISSS = cellfun(@isempty,a_B) ;

if (ISSS(5)==1 && ISSS(6))
else
    error('This method only works when  a_b5 and a_b6 are unprescribed')
end
nn = [1:4] ;
bb = [5:6] ; 
UbB = [R(:,bb)] ;

UbBBAR = M*UbB;
% Matrix Q
coeff = (UbBBAR'*UbB)\UbBBAR';
Q = UbB*coeff ;
IbB = speye(size(Q))-Q ;

%% Degrees of freedom p and s
[~,pB]=licols([UbBBAR ]') ; %
sB = setdiff(1:length(DOFB),pB) ;


ZEROS_M = zeros(size(R,2),length(DOFB)) ;

A = [ IbA  -IbB
    (M*R)' ZEROS_M] ;
 a_Bn  = cell2mat(a_B(nn)) ; 
b = [ R*a_A-R(:,nn)*a_Bn ; 
    (R'*M*R)*a_A] ;

% 
% Choose a set of linearly independent equations
       A_aug = [A,b] ;
      [~,INDEX] =  licols(A_aug') ;
      A = A(INDEX,:) ;
      b = b(INDEX) ;
%
% Now select length(INDEX) linearly ind. columns from A
[~,DOFr] =  licols(A) ;
DOFr = DOFr(:) ;

DOFm = setdiff(1:size(A,2),DOFr) ;
DOFm = DOFm(:) ;
% Ar*dR + Am*dM = b --> dR = inv(Ar)*(b-Am*dM)
G = -A(:,DOFr)\(A(:,DOFm)) ;
G = sparse(G) ;
uBAR = A(:,DOFr)\b ;

DOFAB = [DOFA; DOFB] ;
DOFr = DOFAB(DOFr) ;
DOFm = DOFAB(DOFm) ;



