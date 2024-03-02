%% ACADEMIC TEST 4 %%

x0               = 0.2;
cost.cF          = @(x) (x-1)^2;
cost.gF          = @(x) 2*(x - 1);
constraint.cF{1} = @(x) -x;
constraint.gF{1} = @(x) -1;
s.type           = "fmincon";                           
s.ub             = inf;                                 
s.lb             = -inf;                                 
s.constraintCase = {'INEQUALITY'};
s.maxIter        = 100;