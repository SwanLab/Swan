function  der_Fint_q = EvaluateResidualDerivative(DATA,OPERFE,qR,VAR,qL,MATPRO,VARint_n)

FINT = VAR.FINT ;
qL_plus =   qL + DATA.IntervalEvaluationDerivative ;
tauNON_q_plus= OPERFE.tauNON(qL_plus);
tauNONder_q_plus = OPERFE.tauNONder(qL_plus);

VAR.DISP = [tauNON_q_plus;qR] ;

tauNONallDER_q_RR = speye(length(OPERFE.DOFr),length(OPERFE.DOFr)) ;
tauNONallDER_q_LR = sparse(size(tauNONder_q_plus,1),length(OPERFE.DOFr));
tauNONallDER_q_RL =  sparse(length(OPERFE.DOFr),size(tauNONder_q_plus,2));
tauNONallDER_q = [tauNONder_q_plus,tauNONallDER_q_LR
    tauNONallDER_q_RL, tauNONallDER_q_RR ] ;    % The function below is for evaluating internal forces and residual
VAR.FEXT = VAR.FEXT_extended ; 
[VAR,celastST,FgradST,detFgrad,fINTredNON_mst_DOFl] =  ResidualFromDisplacementsVARnonECM(OPERFE,VAR,MATPRO,DATA,VARint_n,tauNONallDER_q) ;


Fint_plus = VAR.FINT ;

der_Fint_q = (Fint_plus(1)-FINT(1))/DATA.IntervalEvaluationDerivative ;