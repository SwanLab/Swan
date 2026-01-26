function [BstQ,KcLINq,dCqLOC,D_QrotALL] = UpdateMatricesRotation(OPERFE,Qrot) 
% Update rotation-dependent EIFEM operators
% JAHO, 29-OCT-2024, HONEST GREENS, PEDRALBES CENTER, BARCELONA 
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/02_PURE_ROTATION.mlx
% ------------------------------------------------------------------------------------------------------------------
% Diagonal matrix rotation matrices
D_QrotALL = QrotINIall_GLOBAL(Qrot,OPERFE.INDEXsparseROTmat) ;
% -----------------------------------
% Configuration update 
% \BstQ \defeq \DiagC{\BmatIst}  \DiagC{\QrotALL^T}  \LboolCall; 
BstQ = (OPERFE.D_BmatIst*D_QrotALL')*OPERFE.LboolCall ; 
% \KcLINq \defeq \LboolCall^T  \DiagC{\QrotALL}   \DiagC{\KcLIN}    \DiagC{\QrotALL^T}   \LboolCall
KcLINq =  (OPERFE.D_KcLINloc* D_QrotALL')*OPERFE.LboolCall ;
KcLINq = OPERFE.LboolCall'*(D_QrotALL*KcLINq) ; 
% ----------------------------------------------
% COMPUTE ROTATIONAL DISPLACEMENTS (LOCAL)
%\dCqLOC =    \DiagC{\QrotINIall} \XcALL - \DiagC{\QrotALL}  \XcALL
dCqLOC = QrotINIall_GLOBAL(OPERFE.QrotINI,OPERFE.INDEXsparseROTmat)'*OPERFE.XcALL;  
dCqLOC = dCqLOC -D_QrotALL'*OPERFE.XcALL;    
