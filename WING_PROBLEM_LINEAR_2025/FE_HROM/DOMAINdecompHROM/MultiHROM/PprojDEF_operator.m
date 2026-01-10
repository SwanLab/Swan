function  [udef,SS] =  PprojDEF_operator(Vrb,Mintf,u,DATALOC) ;
if nargin == 0
    load('tmp3.mat')
elseif nargin == 3
    DATALOC = [] ; 
end


% PROJECTION ONTO THE SPAN OF rigid body  MODES
% \begin{equation}
% \label{eq:1_3}
%  \SprojRB{e}{} \defeq  \Vrb{e}{} (\Vrb{e}{}^T \Mintf{}{I} \Vrb{}{I})^{-1} \Vrb{}{I}^T \Mintf{}{I}, \hspace{1cm} e=1,2 \ldots \nDOM.
% \end{equation}
% (orthogonality is defined in terms of the geometric mass matrix $\Mintf{}{I}$, as demanded by condition \ref{eq:orth1}). Likewise, the operator defined by
% \begin{equation}
% \label{eq:32qwerw}
%   \SprojDEF{}{I} \defeq \ident -  \SprojRB{}{I}, \hspace{1cm} e=1,2 \ldots \nDOM
% \end{equation}

% coeff =    \Vrb{}{I}^T\Mintf{}{I}  u
coeff = Vrb'*Mintf*u ;
% ScaleMatrix =  (\Vrb{}{I}^T \Mintf{}{I} \Vrb{}{I})
ScaleMatrix =   Vrb'*Mintf*Vrb ;

%  coeff_scaled =  (\Vrb{}{I}^T \Mintf{}{I} \Vrb{}{I})^{-1} \Vrb{}{I}^T \Mintf{}{I} u
coeff_scaled  = ScaleMatrix\coeff ;
udef = u - Vrb*coeff_scaled ;
% Checking whether  udef is zero
normTOTAL = norm(u'*Mintf*u,'fro') ;
TOLrel = 1e-6 ;
DATALOC = DefaultField(DATALOC,'TOLrelSVD',TOLrel) ; % = 
TOLrel = DATALOC.TOLrelSVD ; 
DATALOC.TOL = TOLrel*normTOTAL ; 
DATALOC.RELATIVE_SVD = 0; 
[udef,SS,VV]  =  WSVDT(udef,Mintf,DATALOC)  ;

%  else
%      % Filtering out noisy modes
%      % --------------------------
%      TOLrel = 1e-6 ;
%      TOLabsolute = normTOTAL*TOLrel ;
%    [UU,SS,VV]  =  SVDT(udef,TOLabsolute)  ;
%
%  end



