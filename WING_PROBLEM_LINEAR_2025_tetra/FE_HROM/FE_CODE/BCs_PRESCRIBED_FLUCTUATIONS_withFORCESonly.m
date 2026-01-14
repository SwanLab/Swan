function [G,uBAR,DOFs,DOFm] = BCs_PRESCRIBED_FLUCTUATIONS_withFORCESonly(DOFa,DOFb,R,U_b,gammaB,Rbar,M)
%See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/105_BackPeriodic/05_ONLY_FORCES.mlx
if nargin == 0
    load('tmp.mat')
end 

%%%%%
%   \begin{equation}
%     \Q \defeq \ident -     \R (\Rbar^T \R )^{-1} \Rbar^T   + (\gamma_{b}-1) \U_b  (\U_b^T \Mbar  \U_b)^{-1} \U_b^T \Mbar
%    \end{equation}
% % We begin by writing 
%   \begin{equation}
%     \label{eq:AS*OLL---}
%     \d_{A(\alpha)} =  \Q_{\alpha \alpha}  \d_{B(\alpha)}  +  \Q_{\alpha \beta}  \d_{B(\beta)} 
%   \end{equation}
%  where $\alphaB \subset \{1,2 \ldots M\}$ are the indexes of the linearly independent constraints, and $\betaB$ its complementary set.  These latter indexes can be determined by selecting $7$ rows of $\V = [\R,\U_b]$ such that $\V(\betaB,:)$ is invertible. 
%  \item Since $\Q_{\alpha \alpha}$ is invertible as well, we can solve \refeq{eq:ASFOLL---} for $\d_{B(\alpha)}$: 
% 
% 
% 
nrows = length(DOFa) ; 
ncols =  length(DOFa)  ; 
coeffR = (Rbar'*R)\Rbar' ;
Q = speye(nrows,ncols) -R*coeffR ;  
if ~isempty(U_b)
    coeffU  = (U_b'*M*U_b)\(U_b'*M); 
    Q = Q +(gammaB-1)*U_b*coeffU ; 
    [~,betaB]=licols([R,U_b]') ; %
else
    [~,betaB]=licols(R') ; %
end

alphaB = setdiff(1:length(DOFa),betaB) ;

G = sparse(nrows,ncols) ; 

% 
%  Thus, by setting 
%  
%  \begin{equation}
%   \d_S = \coldos{ \d_{A(\beta)}}{\d_{B(\alpha)}}
%  \end{equation}
% 
%  
%  and 
%  
%  \begin{equation}
%   \d_{M} = \coldos{\d_{A(\alpha)}}{\d_{B(\beta)}}
%  \end{equation}

DOFs = [DOFa(betaB);DOFb(alphaB)]; 
DOFm = [DOFa(alphaB); DOFb(betaB)] ; 


%  \begin{equation}
%   \G \defeq    \matcdos{\T_{\beta \alpha}  \T_{\alpha \alpha}^{-1}  }
% { \T_{\beta \beta} - \T_{\beta \alpha}  \T_{\alpha \alpha}^{-1}    \T_{\alpha \beta}  }
% {\T_{\alpha \alpha}^{-1}}
% {-\T_{\alpha \alpha}^{-1}  \T_{\alpha \beta}} 
% %  \end{equation}
% 


USE_CHEAP_invQaa = 1; 

if USE_CHEAP_invQaa == 0
    invQaa = inv(Q(alphaB,alphaB));
else
    %
    %  \item \textbf{Remark}. \refeq{eq:trainings} may be cast in a format $\Q = \ident - \X \Y$ by making
    %   \begin{equation}
    %     \label{eq:Hgdsdd}
    %     \Q \defeq \ident -     \rowdos{\R  }{(1-\gamma_{b}) \U_b}  \coldos{(\Rbar^T \R )^{-1} \Rbar^T }{(\U_b^T \Mbar  \U_b)^{-1} \U_b^T \Mbar}
    %    \end{equation}%
    % In doing so, we can compute $\Q_{\alpha \alpha}^{-1}$ by using  formula
    %  \begin{equation}
    %  (\ident - \X(\alpha,:) \Y(:,\alpha))^{-1} = \ident +  \X(\alphaB,:) (\ident - \Y(:,\alphaB)   \X(:,\alphaB))^{-1} \Y(:,\alphaB)
    % \end{equation}
    
    if ~isempty(U_b)
        
        X = [R,(1-gammaB)*U_b] ; 
        Y = [coeffR;coeffU] ; 
        
    else
        X = R ;
        Y = coeffR ;
        
    end
    %(\ident - \X(\alpha,:) \Y(:,\alpha))^{-1} = \ident +  \X(\alphaB,:) (\ident - \Y(:,\alphaB)   \X(\alphaB,:))^{-1} \Y(:,\alphaB)
    nx = size(X,2) ; 
     coeffXY = (eye(nx) - Y(:,alphaB)*X(alphaB,:))\Y(:,alphaB) ; 
    invQaa = speye(length(alphaB)) + X(alphaB,:)*coeffXY ;
    
end

%
%G_11 = \T_{\beta \alpha}  \T_{\alpha \alpha}^{-1} 
%
iniROW = 1; 
ifinROW = length(betaB) ; 
iniCOL = 1; 
finCOL = length(alphaB) ; 
G(iniROW:ifinROW,iniCOL:finCOL) = Q(betaB,alphaB)*invQaa; 
%
%G_12 = \T_{\beta \beta} - \T_{\beta \alpha}  \T_{\alpha \alpha}^{-1}    \T_{\alpha \beta} 
%
iniCOL = finCOL+1; 
finCOL = size(G,2) ; 
G(iniROW:ifinROW,iniCOL:finCOL) = Q(betaB,betaB) - Q(betaB,alphaB)*invQaa*Q(alphaB,betaB) ;
% G_21 = \T_{\alpha \alpha}^{-1}
iniROW = ifinROW+1; 
ifinROW = size(G,1);
iniCOL = 1; 
finCOL = length(alphaB) ; 
G(iniROW:ifinROW,iniCOL:finCOL) = invQaa; 

% G_22 = {-\T_{\alpha \alpha}^{-1}  \T_{\alpha \beta}

iniCOL = finCOL+1; 
finCOL = size(G,2) ;  
G(iniROW:ifinROW,iniCOL:finCOL) = -invQaa*Q(alphaB,betaB); 







uBAR = zeros(size(G,1),1) ; 

 