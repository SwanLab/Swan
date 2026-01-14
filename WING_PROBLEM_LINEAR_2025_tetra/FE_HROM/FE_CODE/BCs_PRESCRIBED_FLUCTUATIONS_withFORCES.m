function [G,uBAR,DOFs,DOFm] = BCs_PRESCRIBED_FLUCTUATIONS_withFORCES(DOFa,DOFb,R,U_j,U_r)
%See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/105_BackPeriodic/03_IMPOSED_FORCES.mlx
if nargin == 0
    load('tmp.mat')
end 

% \begin{equation}
%  \V \defeq    \begin{bmatrix}
%          \zero & \U_j & \U_r \\
%          \R   & \U_j  &  \zero
%         \end{bmatrix} 
% \end{equation}
nrows = 2*size(R,1) ; 
ncols = size(R,2) + size(U_j,2) + size(U_r,2) ; 
V = sparse(nrows,ncols) ; 

% V_21 = R   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
iniROW = size(U_j,1)+1; 
finROW = size(V,1) ;

iniCOL = 1;
finCOL = size(R,2) ; 
V(iniROW:finROW,iniCOL:finCOL) = R; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% V_12 = U_j
iniCOL = finCOL+1; 
finCOL  = iniCOL + size(U_j,2)-1 ;

iniROW = 1;
finROW = size(R,1) ; 
V(iniROW:finROW,iniCOL:finCOL) = U_j; 

% V_22 = U_j
iniROW = size(U_j,1)+1; 
finROW = size(V,1) ; 
V(iniROW:finROW,iniCOL:finCOL) = U_j; 

% V_13  = U_r
iniROW = 1;
finROW = size(R,1) ; 

iniCOL = finCOL+1; 
finCOL  = size(V,2) ;
V(iniROW:finROW,iniCOL:finCOL) = U_r; 
  
 
 

% 
 %\item Let $\alphaB$ and $\betaB$ be two sets of rows of $\V$ such that $\V_{\alpha} = \V(\alphaB,:)$ is invertible 
 %and $\betaB \union \alphaB = 1,2 \ldots 2 M$.
 % We shall asign the role of master  and slave DOFs to the entries of $\d_F$ corresponding to $\alphaB$ and $\betaB$, respectively:  
%  \begin{equation}
%   \d_M  = \d_{F(\alpha)}
%  \end{equation}
%  \begin{equation}
%   \d_S = \d_{F(\beta)}
%  \end{equation}
 
[~,alphaB]=licols(V') ; %
 
betaB = setdiff(1:size(V,1),alphaB) ;
betaB = betaB(:) ;
DOFf = [DOFa;DOFb]; 
DOFs = DOFf(betaB) ; 
DOFm = DOFf(alphaB) ; 

% \begin{equation}
%  \G \defeq \V_{\beta} \V_{\alpha}^{-1}
% \end{equation}
G = V(betaB,:)*inv(V(alphaB,:)) ; 

 


uBAR = zeros(size(G,1),1); 



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

 