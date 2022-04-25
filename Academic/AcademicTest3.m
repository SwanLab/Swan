%% ACADEMIC TEST 3 %%

x0                           = [3,3];
j.cF                         = @(x) x(1)^2 + (x(2)+3)^2;                 % COST
j.gF                         = @(x) [2*x(1); 2*(x(2) + 3)];              % GRADIENT COST
c.cF                         = @(x) [-x(1)^2 + x(2),-x(1)-x(2)-2];       % CONST
c.gF                         = @(x) [-2*x(1), -1; 1, -1];                % GRADIENT CONST
s.type                       = "fmincon";                                % CONST OPTIMIZER
s.uncOptimizerSettings.ub    = inf;                                      % UPPER BOUND
s.uncOptimizerSettings.lb    = -inf;                                     % LOWER BOUND