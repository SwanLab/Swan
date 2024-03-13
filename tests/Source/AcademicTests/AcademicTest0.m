%% ACADEMIC TEST 0 %%

x0               = [3;2];
cost.cF          = @(x) x(1)^2 + 2 * x(2)^2;
cost.gF          = @(x) [2*x(1);4*x(2)];
constraint.cF{1} = @(x) 2*x(1) + x(2) - 9;
constraint.gF{1} = @(x) [2; 1];
constraint.cF{2} = @(x)  x(1) + 2*x(2) - 10;
constraint.gF{2} = @(x) [1; 2];
s.type           = "fmincon";                                
s.ub             = inf;                                    
s.lb             = -inf;                                   
s.constraintCase = {'INEQUALITY','EQUALITY'};
s.maxIter        = 100;