function [BstQ,KcLINq,D_QrotALL] = UpdateMatricesRotationV2(OPERFE,Qrot) 
% Update rotation-dependent EIFEM operators
% JAHO, 02-DEC-2024, HONEST GREENS, PEDRALBES CENTER, BARCELONA 
% Latex notation:: /home/joaquin/Desktop/CURRENT_TASKS/PAPERS_2020_onwards/12_EIFEM_EXTENS/EIFEM_largeROTexAPPV.tex 
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/04_UPDrotUNC.mlx
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


% % COMPUTE ROTATIONAL DISPLACEMENTS (LOCAL)
% %\dCqLOC =    \DiagC{\QrotINIall} \XcALL - \DiagC{\QrotALL}  \XcALL
% dCqLOC = QrotINIall_GLOBAL(OPERFE.QrotINI,OPERFE.INDEXsparseROTmat)'*OPERFE.XcALL;  
% dCqLOC = dCqLOC -D_QrotALL'*OPERFE.XcALL;    
