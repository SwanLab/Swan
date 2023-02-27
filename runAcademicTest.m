clear
close all

filename = 'AcademicTest1';
run(filename);
d                            = DesignVariableAcademic();
d.init(x0);
j.dV                         = d;
c.dV                         = d;
s.designVar                  = d;
s.cost                       = AcademicCost(j);
s.constraint                 = AcademicConstraint(c);
s.constraint.nSF             = nConstr;
s.nConstraints               = nConstr;
s.dualVariable               = DualVariable(s);
s.outputFunction.type        = "Academic";
s.outputFunction.iterDisplay = "iter";
s.outputFunction.monitoring  = MonitoringManager(s);
s.optimizerNames.primal     = 'PROJECTED GRADIENT';
opt = Optimizer.create(s);
opt.solveProblem();
result = d.value;