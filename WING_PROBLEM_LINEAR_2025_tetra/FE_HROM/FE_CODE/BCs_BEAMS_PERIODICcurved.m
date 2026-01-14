function [G,uBAR,DOFr,DOFm] = BCs_BEAMS_PERIODICcurved(DOFA,DOFB,R,a_A,a_B,Rbar,Qsingle)

if nargin == 0
    load('tmp1.mat')
end



nconstr = length(a_A) + length(DOFA) ;

% Assembly of Q diagonal 
% ----------------------
%
ndim = size(Qsingle,1) ; 
nnode = size(R,1)/ndim ; 
Q =cell(nnode,1) ; 
Q(:) = {sparse(Qsingle)} ; 
Q = blkdiag(Q{:}) ;
%%  \S = \Q \Mbar \R  = \Q \Rbar , where Q = diag(Q,Q,Q,Q)
S = Q*Rbar ; 
%Sother = RotateMatrix(Q,Rbar) ;  % Other possible method


% Select 6 linearly independent rows from S
[~,r]=licols(S') ; %
l = setdiff(1:length(DOFA),r) ;
% 
da = a_A-a_B ; 
%% Computation of J  and b
% \J \defeq -\S_r^{-T} \S_l^T
% \b \defeq \S_r^{-T}    \Rbar^T\R \a_B
J = -S(r,:)'\S(l,:)' ; 
b = S(r,:)'\(Rbar'*(R*a_B)) ;
%%- --------------------------------------
%  \H = \defeq \Q_r^T \J + \Q_l^T
%  \c \defeq \Q_r^T \b + \R \da
H = Q(r,:)'*J + Q(l,:)' ; 
c = Q(r,:)'*b + R*da ; 

%
% \d_R = \G \d_M + \uBAR \RRaam 
% \coldos{\d_A}{\d_{Br}} = \coldos{\H}{\J} \d_{Bl} + \coldos{\c}{\b} 


DOFr = [DOFA; DOFB(r)  ]  ; 
DOFm = DOFB(l) ; 
G = sparse(length(DOFr),length(DOFm)) ; 
uBAR = zeros(length(DOFr),1) ; 

iini = 1; 
ifin = size(H,1) ; 
G(iini:ifin,:) = H ; 
uBAR(iini:ifin) = c  ; 

iini = ifin+1; 
ifin = iini+size(J,1)-1 ; 
G(iini:ifin,:) = J ; 
uBAR(iini:ifin) = b ; 
 




 


%  
% 
% 
% nrows = length(indx) + 3 ;
% ncols = length(indyz) -3  ;  %
% G_A = sparse(nrows, ncols) ;
% G_B = sparse(nrows, ncols) ;
% uA = zeros(nrows,1) ;
% uB = zeros(nrows,1) ;
% 
% 
% Rx = R(indx,:) ;
% COLn = [2,3,4] ;
% Ryz = R(indyz,:) ;
% Ryz_n = R(indyz,COLn) ;
% % Select 3 linearly independent rows from Ryz_n
% [~,r]=licols(Ryz_n') ; %
% l = setdiff(1:length(indyz),r) ;
% R_r = Ryz_n(r,:) ;
% R_l = Ryz_n(l,:) ;
% DOFr_A = [DOFA(indx); DOFA(indyz(r))] ; % Slave DOFS
% DOFr_B = [DOFB(indx); DOFB(indyz(r))] ;
% DOFm_A = DOFA(indyz(l)) ; % Master DOFs
% DOFm_B = DOFB(indyz(l)) ;
% %
% ROWS =  (length(indx)+1):(length(indx)+3) ;
% G_A(ROWS,:) = -inv(R_r')*R_l' ;
% G_B(ROWS,:) = -inv(R_r')*R_l' ;
% %
% uA(1:length(indx)) = Rx*a_A ;
% uB(1:length(indx)) = Rx*a_B ;
% 
% uA(length(indx)+1:end) = inv(R_r')*(Ryz_n'*(Ryz*a_A));
% uB(length(indx)+1:end) = inv(R_r')*(Ryz_n'*(Ryz*a_B));
% 
% DOFr = [DOFr_A; DOFr_B] ;
% DOFm = [DOFm_A; DOFm_B] ;
% 
% G = sparse(length(DOFr),length(DOFm) ) ;
% 
% ROWS = 1:length(DOFr_A) ; COLS  = 1:length(DOFm_A)  ;
% G(ROWS,COLS) = G_A ;
% ROWS = (length(DOFr_A)+1 ):length(DOFr); COLS  = (length(DOFm_A)+1):length(DOFm)  ;
% G(ROWS,COLS) = G_B ;
% 
% 
% uBAR = [uA; uB] ;
% 
% % Sorting for avoiding large bandwidth
% [DOFr,INDR ]= sort(DOFr) ;
% [DOFm,INDM ]= sort(DOFm) ;
% G = G(INDR,INDM) ;
% uBAR = uBAR(INDR) ;
% 
% 
% 
