%% ACADEMIC TEST 3 %%

x0                           = [3;3];
j.cF                         = @(x) (x(1)+3)^2 + (x(2))^2;                 % COST
j.gF                         = @(x) [2*(x(1)+3); 2*(x(2))];              % GRADIENT COST
c.cF                         = @(x) [-x(1)^2 + x(2);-x(1)-x(2)-2];       % CONST
c.gF                         = @(x) [-2*x(1), -1; 1, -1]';                % GRADIENT CONST
s.type                       = "IPM";                                % CONST OPTIMIZER
s.uncOptimizerSettings.ub    = inf;                                      % UPPER BOUND
s.uncOptimizerSettings.lb    = -inf;                                     % LOWER BOUND
nConstr                      = 2;                                        % NUMBER OF CONSTRAINTS
s.incrementalScheme.iStep    = 1;
s.incrementalScheme.nSteps   = 1;
s.dualVariable               = [];
s.maxIter                    = [];
s.targetParameters.constr_tol = 1e-3;
s.constraintCase             = {'INEQUALITY','INEQUALITY'};
s.maxIter                    = 1000;
s.postProcessSettings.shallPrint = false;