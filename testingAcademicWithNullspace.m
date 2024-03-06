close all;
clear;

%% Problem 1

% Min problem
cost.cF = @(x) x(1) + x(2) + x(3);
cost.gF = @(x) [1; 1; 1];

constraint.cF{1} = @(x) x(1).^2 + x(2).^2-2;
constraint.gF{1} = @(x) [2*x(1); 2*x(2); 0];

constraint.cF{2} = @(x) x(1)+x(3)-1;
constraint.gF{2} = @(x) [1; 0; 1];

% Solution
% x1=0; x2=-sqrt(2)=-1.4142; x3=1; lambda1=0.25*sqrt(2)=0.3536; lambda2=-1


% Setting up
x0               = [5;5;5];
s.type           = "NullSpace";
s.ub             = inf;
s.lb             = -inf;
s.maxIter        = 400;
s.constraintCase = {'EQUALITY','EQUALITY'};
s.etaNorm        = 0.02;
s.gJFlowRatio    = 1;

cParams.cost         = cost;
cParams.constraint   = constraint;
cParams.initialGuess = x0;
cParams.settings     = s;
problem              = AcademicProblem(cParams);

problem.compute();
xStar = problem.result;