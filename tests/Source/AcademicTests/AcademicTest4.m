%% ACADEMIC TEST 4 %%

x0 = 0.2;
cost.cF = @(x) (x-1)^2;
cost.gF = @(x) 2*(x - 1);
constraint.cF = @(x) -x;
constraint.gF = @(x) -1;
s.type                       = "fmincon";                                % CONST OPTIMIZER
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
s.shallPrint = false;