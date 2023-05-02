%% ACADEMIC TEST 5 %%

x0 = 0.2;
j.cF = @(x) (x(1)-1)^2;
j.gF = @(x) 2*(x(1) - 1);
j.hF = @(x) 2;
c.cF = @(x) -x;
c.gF = @(x) -1;
s.type                       = "IPM";                                % CONST OPTIMIZER
s.uncOptimizerSettings.ub    = inf;                                      % UPPER BOUND
s.uncOptimizerSettings.lb    = -inf;                                     % LOWER BOUND
nConstr                      = 1;                                        % NUMBER OF CONSTRAINTS
s.incrementalScheme.iStep    = 1;
s.incrementalScheme.nSteps   = 1;
s.dualVariable               = [];
s.maxIter                    = [];
s.targetParameters.constr_tol = 1e-3;
s.constraintCase             = {'INEQUALITY'};
s.maxIter                    = 100;
s.postProcessSettings.shallPrint = false;