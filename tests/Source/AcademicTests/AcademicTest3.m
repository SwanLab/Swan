%% ACADEMIC TEST 3 %%

x0                           = [3;3];
cost.cF                         = @(x) (x(1)+3)^2 + (x(2))^2;                 % COST
cost.gF                         = @(x) [2*(x(1)+3); 2*(x(2))];              % GRADIENT COST
constraint.cF                         = @(x) [-x(1)^2 + x(2);-x(1)-x(2)-2];       % CONST
constraint.gF                         = @(x) [-2*x(1), -1; 1, -1]; 
s.type                       = "fmincon";                                % CONST OPTIMIZER
s.uncOptimizerSettings.ub    = inf;                                      % UPPER BOUND
s.uncOptimizerSettings.lb    = -inf;                                     % LOWER BOUND
nConstr                      = 2;                                        % NUMBER OF CONSTRAINTS
s.incrementalScheme.iStep    = 1;
s.incrementalScheme.nSteps   = 1;
s.dualVariable               = [];
s.maxIter                    = [];
s.targetParameters.constr_tol = 1e-3;
s.constraintCase             = {'INEQUALITY','INEQUALITY'};
s.maxIter                    = 100;
s.postProcessSettings.shallPrint = false;
s.shallPrint = true;