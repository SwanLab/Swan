function [pST,kPRESSst,Y_neg]=  HydrostaticPressure(OPERHYDRO,VAR,DATA) ;
% Hydrostatic pressure (waterline X2 = 0)
% % See formulation in /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/DEVELOPMENTS/DynamicFloatingBodies.tex
%  JAHO, 1-July-2021, Cartagena, SPAIN (UPCT)
if nargin == 0
    load('tmp.mat')
end

% % \begin{equation}
%   \boxed{\pPRESSst =  -g \densW \HEAV{-\yPRESSst} \hadmP \yPRESSst   }
%  \end{equation}

g = -DATA.vGRAVITY(2) ;
densW = DATA.FOLLOWER_LOADS.HYDROSTATIC.densW ;


%  \yPRESSst & = \yPRESSiniST +  \diag{({\NbSTgrav{1}  },{\NbSTgrav{2}  },{\cdots},{\NbSTgrav{\nelemB}} )} \LLbB{}{} \d

% A) Determination of the x2 coordinate of all the Gauss points of the
% boundary
% --------------------------------------------------------------------
dyPRESSst=OPERHYDRO.yPRESSiniST +  ConvertBlockDiag_general(OPERHYDRO.NbSTgrav,1,OPERHYDRO.irows1,OPERHYDRO.icols1)*OPERHYDRO.Lbool*VAR.DISP(:) ;
% Y_neg= \HEAV{-\yPRESSst} \hadmP \yPRESSst
Y_neg = dyPRESSst.*(dyPRESSst<0);

pST = -g*densW*Y_neg ;

% Contribution to the stiffness matrix
% ----------------------------------------
% \begin{equation}
% \label{eq:73dgwedd}
%   \kPRESSst^e    \defeq  \colcuatro{{\kPRESS^e}_{\xi_1}}{{\kPRESS^e}_{\xi_2}}{\vdots}{{\kPRESS^e}_{\xi_k}}  = g \densW \HEAV{-\yPRESSst^e} \hadmP   \NbSTgrav{e}
% \end{equation}

if    DATA.TMP.NO_COMPUTE_STIFFNESS_HYDRO == 0
    kPRESSst = g*densW*bsxfun(@times,OPERHYDRO.NbSTgrav,(dyPRESSst<0)) ;
    
else
    kPRESSst = 0 ;
end
