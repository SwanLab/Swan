function  [AsnapDISPdef,AsnapDISPrb,qRB] =  SprojDEF_operator(PhiRB,Mdom,AsnapDISP) ; 
if nargin == 0
    load('tmp.mat')
end

 
% PROJECTION ONTO THE SPAN OF rigid body  MODES   
% \begin{equation}
% \label{eq:1_3}
%  \SprojRB{e}{} \defeq  \PhiRB{e}{} (\PhiRB{e}{}^T \Mdom{e}{} \PhiRB{e}{})^{-1} \PhiRB{e}{}^T \Mdom{e}{}, \hspace{1cm} e=1,2 \ldots \nDOM. 
% \end{equation}
% (orthogonality is defined in terms of the geometric mass matrix $\Mdom{e}{}$, as demanded by condition \ref{eq:orth1}). Likewise, the operator defined by  
% \begin{equation}
% \label{eq:32qwerw}
%   \SprojDEF{e}{} \defeq \ident -  \SprojRB{e}{}, \hspace{1cm} e=1,2 \ldots \nDOM 
% \end{equation} 

% coeff =    \PhiRB{e}{}^T \Mdom{e}{}  AsnapDISP  

if isempty(Mdom)
    Mdom = speye(size(PhiRB,1)) ; 
end


coeff = PhiRB'*Mdom*AsnapDISP ; 
% ScaleMatrix =  (\PhiRB{e}{}^T \Mdom{e}{} \PhiRB{e}{}) 
ScaleMatrix =   PhiRB'*Mdom*PhiRB ; 

%  coeff_scaled =  (\PhiRB{e}{}^T \Mdom{e}{} \PhiRB{e}{})^{-1} \PhiRB{e}{}^T \Mdom{e}{} AsnapDISP
qRB  = ScaleMatrix\coeff ; 

% Self-equilibrated part
AsnapDISPrb =   PhiRB*qRB ; 

AsnapDISPdef = AsnapDISP -AsnapDISPrb ; 

% Checking whether  AsnapDISPdef is zero 
normTOTAL = norm(AsnapDISP,'fro') ; 
normDEF = norm(AsnapDISPdef,'fro') ; 


 rel_norm = normDEF/normTOTAL ; 
 TOL = 1e-12 ; 
 if rel_norm <=TOL
     AsnapDISPdef = [] ; % The projection is zero
 end



