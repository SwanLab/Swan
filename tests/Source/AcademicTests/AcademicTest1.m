%% ACADEMIC TEST 1 %%

x0                           = [1.5;2.25];
cost.cF                         = @(x) 0.3*x(1) + x(2);                     % COST
cost.gF                         = @(x) [0.3; 1];                            % GRADIENT COST
constraint.cF                         = @(x) [1/x(1) - x(2);x(1)+x(2)-3];         % CONST
constraint.gF                         = @(x) [-1/x(1)^2, 1; -1, 1];               % GRADIENT CONST
s.type                       = "fmincon";                                % CONST OPTIMIZER
s.uncOptimizerSettings.ub    = inf;                                      % UPPER BOUND
s.uncOptimizerSettings.lb    = -inf;                                     % LOWER BOUND
nConstr                      = 2;                                        % NUMBER OF CONSTRAINTS
s.incrementalScheme.iStep    = 1;
s.incrementalScheme.nSteps   = 1;
s.dualVariable               = [];
s.maxIter                    = [];
s.targetParameters           = [];
s.constraintCase             = {'INEQUALITY','INEQUALITY'};
s.maxIter                    = 50;
s.postProcessSettings.shallPrint = false;
s.shallPrint = true;