%% ACADEMIC TEST 2 %%

x0                           = [1.5,2.25];
j.cF                         = @(x) (x(1)-2)^2 + (x(2)-2)^2;             % COST
j.gF                         = @(x) [0.3; 1];                            % GRADIENT COST
c.cF                         = @(x) [1/x(1) - x(2),x(1)+x(2)-3];         % CONST
c.gF                         = @(x) [-1/x(1)^2, 1; -1, 1];               % GRADIENT CONST
s.type                       = "fmincon";                                % CONST OPTIMIZER
s.uncOptimizerSettings.ub    = inf;                                      % UPPER BOUND
s.uncOptimizerSettings.lb    = -inf;                                     % LOWER BOUND