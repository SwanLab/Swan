function [Q_d,Q_r] = ObjectiveFunction_INTF(qRB,Pcomp,INDrig,INDdef,qDEF,bCOMP,mCOMP,Treac,rDEF)

Q_d = (qRB'*Pcomp(INDrig,INDrig)*qRB + 2*qRB'*Pcomp(INDrig,INDdef)*qDEF+qDEF'*Pcomp(INDdef,INDdef)*qDEF...
    -2*qRB'*bCOMP(INDrig)-2*qDEF'*bCOMP(INDdef) + ...
    mCOMP) ;
Q_r = (rDEF'*Treac*rDEF ) ;
disp('--------------------------------------')
disp(['Objective function DISPLACEMENT'])
disp(['Q_q =',num2str((Q_d))])
disp('--------------------------------------')
disp(['Objective function REACTIONS'])
disp(['Q_q =',num2str((Q_r))])