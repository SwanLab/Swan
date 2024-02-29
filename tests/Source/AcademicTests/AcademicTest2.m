%% ACADEMIC TEST 2 %%

x0               = [1.5;2.25];
cost.cF          = @(x) (x(1)-2)^2 + (x(2)-2)^2;             % COST
cost.gF          = @(x) [2*x(1) - 4; 2*x(2) - 4];            % GRADIENT COST
constraint.cF{1} = @(x) 1/x(1) - x(2);                       % CONST
constraint.gF{1} = @(x) [-1/x(1)^2; -1];                     % GRADIENT CONST
constraint.cF{2} = @(x) x(1)+x(2)-3;                         % CONST
constraint.gF{2} = @(x) [1; 1];                              % GRADIENT CONST
s.type           = "fmincon";                                % CONST OPTIMIZER
s.ub             = inf;                                      % UPPER BOUND
s.lb             = -inf;                                     % LOWER BOUND
s.constraintCase = {'INEQUALITY','INEQUALITY'};
s.maxIter        = 50;