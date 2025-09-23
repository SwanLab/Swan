function LHS = IntegrateLHS(f,test,trial,mesh,quadOrder)
if nargin < 5 || isempty(quadOrder)
    qTe = test.getOrderNum();
    qTr = trial.getOrderNum();
    quadOrder = qTe + qTr;
end
s.trial = trial;
s.test  = test;
s.mesh  = mesh;
s.quadratureOrder = quadOrder;
lhs = LHSIntegrator(s);
LHS = lhs.compute(f);
end