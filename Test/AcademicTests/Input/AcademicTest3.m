%% ACADEMIC TEST 3 %%

x0               = [3;3];
cost.cF          = @(x) (x(1)+3)^2 + (x(2))^2;               % COST
cost.gF          = @(x) [2*(x(1)+3); 2*(x(2))];              % GRADIENT COST
constraint.cF{1} = @(x) -x(1)^2 + x(2);                      % CONST
constraint.gF{1} = @(x) [-2*x(1); 1]; 
constraint.cF{2} = @(x) -x(1)-x(2)-2;                        % CONST
constraint.gF{2} = @(x) [-1; -1]; 
s.type           = "fmincon";                                % CONST OPTIMIZER
s.ub             = inf;                                      % UPPER BOUND
s.lb             = -inf;                                     % LOWER BOUND
s.constraintCase = {'INEQUALITY','INEQUALITY'};
s.maxIter        = 100;