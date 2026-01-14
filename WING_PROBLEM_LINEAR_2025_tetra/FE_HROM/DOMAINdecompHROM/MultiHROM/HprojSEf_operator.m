function [AsnapREACse,MdomFFinv] = HprojSEf_operator(PsiRBf,MdomFF,AsnapREAC) 
if nargin == 0
    load('tmp3.mat')
end

%     In analogy to \refeq{eq:1_3} and \refeq{eq:32qwerw}, we introduce the orthogonal projection matrices
%  \begin{equation}
%  \label{eq:alter14}
%   \HprojRESf{e}{} \defeq \PsiRBf{e}{} (\PsiRBf{e^T}{} \MdomFF{e^{-1}}{} \PsiRBf{e}{} )^{-1}  {\PsiRBf{e^T}{ } \MdomFF{e^{-1}}{}},  \hspace{1cm}     e= 1,2 \ldots \nDOM
%  \end{equation}
%  and 
%  \begin{equation}
%    \label{eq:alter15}
%      \HprojSEf{e}{}  = \ident - \HprojRESf{e}{},  \hspace{1cm}     e= 1,2 \ldots \nDOM
%  \end{equation}


% PROJECTION ONTO THE SPAN OF RESULTANT MODES PsiRBf 
% Inverse MdomFF 

MdomFFinv = inv(MdomFF) ; 
% MdomFFinv{1,1} = inv(MdomFF{1,1}) ; 
% MdomFFinv{2,2} = inv(MdomFF{2,2}) ;

% coeff = {\PsiRBf{e^T}{ } \MdomFF{e^{-1}}{}} AsnapREAC
coeff =  (PsiRBf)'*(MdomFFinv)*AsnapREAC ; 
% ScaleMatrix = (\PsiRBf{e^T}{} \MdomFF{e^{-1}}{} \PsiRBf{e}{} )^{-1}
ScaleMatrix =  (PsiRBf)'*(MdomFFinv)*(PsiRBf) ; 

%  coeff_scaled = (\PsiRBf{e^T}{} \MdomFF{e^{-1}}{} \PsiRBf{e}{} )^{-1}  {\PsiRBf{e^T}{ } \MdomFF{e^{-1}}{}} AsnapREAC
coeff_scaled  = ScaleMatrix\coeff ; 

% Self-equilibrated part
AsnapREACse = AsnapREAC -  (PsiRBf)*coeff_scaled ; 



