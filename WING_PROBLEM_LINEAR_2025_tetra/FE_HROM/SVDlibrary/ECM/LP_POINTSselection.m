function   [setPoints,WdomRED,fval,exitflag,output,lambda] = LP(SNAPfNW,W,TOLv,TOLg,TOLSVD) ; 
% Selection method from a set of Gauss points using Linear Programmming 
% SNAPfNW --> Matrix of Snapshots 
% W - >Weights 
% TOLv --> Tolerance for integrating the volume
% TOLg -->Tolerance for integrating the rest of the functions 
% See LP_basedmethod.pdf 
% JAHO 

if nargin==0
    load('tmp2.mat')
end

% ``Truth" integral 
INTexac  = SNAPfNW.'*W ; 
V = sum(W) ; 
M = length(W) ; 
Xf = zeros(size(SNAPfNW)) ; 
for imode = 1:size(SNAPfNW,2)
    iii =  (SNAPfNW(:,imode)-INTexac(imode)/V) ;  
    Xf(:,imode) = iii; 
end
% Now Xf is such that Xf.'*W  = 0 
% Now we apply the SVD 
[Gt,~,~] = RSVDT(Xf,TOLSVD) ; 
% U is the matrix containing the constraints. 
% ECM ---> 
G  = Gt' ; % [U,sqrt(W)];
n = size(G,1) ; 


%%% LINEAR PROGRAM 
UNOS = ones(M,1) ; 
f = UNOS;   % Objective function 
% Constraint matrices 

% Number of constraints 
nconst = 2*(1+size(G,1))  + M ; 
A = sparse(nconst,M) ; 
b = zeros(nconst,1) ; 

% Positive weights 
iini = 1; 
ifin = iini + M -1; 
A(iini:ifin,:) = -speye(M) ; 
b(iini:ifin) = 0; 
% Volume const. 
% Pos
iini = ifin +1 ;
ifin = iini ; 
A(iini:ifin,:) = UNOS' ; 
b(iini:ifin)  = (1+TOLv)*V; 
% Neg 
iini = ifin +1 ;
ifin = iini ; 
A(iini:ifin,:) = -UNOS' ; 
b(iini:ifin)  = (-1+TOLv)*V;
% Function const.  
% Pos (n)
iini = ifin +1 ;
ifin = iini + n-1 ; 
A(iini:ifin,:) =  G ; 
b(iini:ifin)  =  TOLg;
% Neg (n)
iini = ifin +1 ;
ifin = iini + n-1 ; 
A(iini:ifin,:) =  -G ; 
b(iini:ifin)  =  TOLg;

%%% 
options = optimoptions('linprog','Algorithm','dual-simplex');
x0 = [];
Aeq = []; beq = []; lb = []; ub = [];
 
 
disp('Invoking linprog ...')
[w,fval,exitflag,output,lambda] = linprog(f,A,b,Aeq,beq,lb,ub,x0,options) ; 
disp('...finished') 
 

[setPoints  bbbb] = find(w>0) ; 

WdomRED = w(setPoints) ; 




 