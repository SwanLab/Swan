function [G,uBAR,DOFs,DOFm] = BCs_PRESCRIBED_FLUCTUATIONS(DOFa,DOFb,R,a_A,a_B,Rbar,Mbar,U)

if nargin == 0
    load('tmp.mat')
end 

disp(['Checking that the columns of  U are mutually  orthogonal '])

Check1  =U'*Mbar*U - eye(size(U,2)) ;
if norm(Check1,'fro') > 1e-10
    error('U is not orthogonal ')
else
    disp(['....done'])
end

disp(['Checking that U is M-orthogonal to the rigid body modes '])

Check2  =U'*Mbar*R  ;
if norm(Check2,'fro') > 1e-10
    error('U is not orthogonal to R  ')
else
    disp(['....done'])
end

%   Let $\alphaB$ and $\betaB$ be a set of rows of $\U$ such that $\U_{\alpha} = \U(\alphaB,:)$ is invertible and $\betaB
% \union \alphaB = \{1,2 \ldots N\}$, 

 [~,alphaB]=licols(U') ; %
betaB = setdiff(1:length(DOFa),alphaB) ;
% 
%   \begin{equation}
%    \W \defeq  \U_{\beta} \U_{\alpha}^{-1}   
%   \end{equation}

W = U(betaB,:)*inv(U(alphaB,:)) ; 

%   \begin{equation}
%     \Y =  \R(\beta,:) - \U_{\beta} \U_{\alpha}^{-1}  \R(\alpha,:)
%   \end{equation}

Y = R(betaB,:) - W* R(alphaB,:) ;  

% Slave DOFs/Master  DOFs
% 
% \begin{equation}
%  \coldos{ \d_{A(\beta)}}{ \d_{B(\beta)}} = \matcdos{\W}{\zero}{\zero}{\W}  \coldos{ \d_{A(\alpha)}}{ \d_{B(\alpha)}}  +  \coldos{\Y \a_A}{\Y \a_B} 
% \end{equation}

DOFs = [DOFa(betaB); DOFb(betaB)] ; 
DOFm = [DOFa(alphaB); DOFb(alphaB)] ; 

G = sparse(length(DOFs),length(DOFm)) ; 
uBAR = zeros(length(DOFs),1) ; 

iniROW = 1; 
finROW = size(W,1) ; 
iniCOL = 1; 
finCOL = size(W,2) ;

G(iniROW:finROW,iniCOL:finCOL) = W ; 
uBAR(iniROW:finROW) = Y*a_A ; 

iniROW = finROW+1; 
finROW = size(G,1) ; 
iniCOL = finCOL+1; 
finCOL = size(G,2) ;

G(iniROW:finROW,iniCOL:finCOL) = W ; 
uBAR(iniROW:finROW) = Y*a_B ; 



% 
%  iini = 1; 
% ifin = size(J,1) ; 
% G(iini:ifin,:) = - J ; 
% uBAR(iini:ifin) = b  ; 
% 
% iini = ifin+1; 
% ifin = iini+size(J,1)-1 ; 
% G(iini:ifin,:) = - J ; 
% uBAR(iini:ifin) = b - R(r,:)*da; 
% 
% iini = ifin +1; 
% G(iini:end,:) = speye(length(l)); 
% uBAR(iini:end) =   - R(l,:)*da; 
% 

 