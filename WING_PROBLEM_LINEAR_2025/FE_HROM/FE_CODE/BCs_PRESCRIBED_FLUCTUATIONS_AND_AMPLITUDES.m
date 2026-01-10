function [G,uBAR,DOFs,DOFm] = BCs_PRESCRIBED_FLUCTUATIONS_AND_AMPLITUDES(DOFa,DOFb,R,a_A,a_B,Rbar,Mbar,U,D)

if nargin == 0
    load('tmp.mat')
end 

% 
%  \begin{equation}
%   \V \defeq \coldos{\U \D^A }{\U \D^B}
%  \end{equation}
V = [U*D{1};U*D{2}] ; 
% \item In general, $\V$ will not be full rank, so it is necessary to replace it by a basis matrix for its column space:
 
%  \begin{equation}
%   [\V,\bullet,\bullet]  \leftarrow \SVDT{\V}
%  \end{equation}
[V,~,~] = SVDT(V) ; 

% 
%  
%  \item Let $\alphaB$ and $\betaB$ be two sets of rows of $\V$ such that $\V_{\alpha} = \V(\alphaB,:)$ 
%is invertible and $\betaB \union \alphaB = 1,2 \ldots 2 M$.  
%We shall asign the role of master DOFs to the entries of $\d_F$ corresponding to $\alphaB$, 
%and the slave DOFs the complementary set:  
%  \begin{equation}
%   \d_M  = \d_{F(\alpha)}
%  \end{equation}
%  \begin{equation}
%   \d_S = \d_{F(\beta)}
%  \end{equation}
if size(V,2) == 1
    [~,alphaB ] = max(abs(V)) ; 
else
[~,alphaB]=licols(V') ; %
end
betaB = setdiff(1:size(V,1),alphaB) ;
betaB = betaB(:) ;
DOFf = [DOFa;DOFb]; 
DOFs = DOFf(betaB) ; 
DOFm = DOFf(alphaB) ; 

% \begin{equation}
%  \G \defeq \V_{\beta} \V_{\alpha}^{-1}
% \end{equation}
G = V(betaB,:)*inv(V(alphaB,:)) ; 

% %  \begin{equation}
%   \y \defeq \coldos{\R \a_A}{ \R \a_B} 
%  \end{equation}
%  


y = [R*a_A; R*a_B] ; 

% \begin{equation}
%  \uBARb \defeq \y_{\beta} -  \V_{\beta}  \V_{\alpha}^{-1} \y_{\alpha}
% \end{equation}


uBAR = y(betaB) - G*y(alphaB) ; 



% G = sparse(length(DOFs),length(DOFm)) ; 
% uBAR = zeros(length(DOFs),1) ; 
% 
% iniROW = 1; 
% finROW = size(W,1) ; 
% iniCOL = 1; 
% finCOL = size(W,2) ;
% 
% G(iniROW:finROW,iniCOL:finCOL) = W ; 
% uBAR(iniROW:finROW) = Y*a_A ; 
% 
% iniROW = finROW+1; 
% finROW = size(G,1) ; 
% iniCOL = finCOL+1; 
% finCOL = size(G,2) ;
% 
% G(iniROW:finROW,iniCOL:finCOL) = W ; 
% uBAR(iniROW:finROW) = Y*a_B ; 



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

 