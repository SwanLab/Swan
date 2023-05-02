%% ACADEMIC TEST 4 %%

x0 = [-7;-5];
j.cF = @(x) x(1)^2-2*x(1)*x(2)+4*x(2)^2;
j.gF = @(x) [2*x(1)-2*x(2);-2*x(1)+8*x(2)];
j.hF = @(x) [2,-2;-2,8];
c.cF = @(x) -[0.1*x(1)-x(2)-1;-10*x(1)+x(2)-1];
c.gF = @(x) -[0.1,-1.0;-10.0,1.0];
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
s.maxIter                    = 100;
s.postProcessSettings.shallPrint = false;