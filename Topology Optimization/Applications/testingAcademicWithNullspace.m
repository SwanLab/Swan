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
cParams.printingPath = true;
problem              = AcademicProblem(cParams);

problem.compute();
xStar = problem.result;

%% Problem 2

close all;
clear;

% Min problem
cost.cF = @(x) exp(x(1).*x(2));
cost.gF = @(x) [x(2).*exp(x(1).*x(2)); x(1).*exp(x(1).*x(2))];

constraint.cF{1} = @(x) x(1).^2 + x(2).^2-8;
constraint.gF{1} = @(x) [2*x(1); 2*x(2)];

% Solutions
% x1=2; x2=-2;
% x1=-2; x2=2;


% Setting up
x0               = [1;0];
s.type           = "NullSpace";
s.ub             = inf;
s.lb             = -inf;
s.maxIter        = 2000; % - - 206 inf
s.constraintCase = {'EQUALITY'};
s.etaNorm        = 0.01;
s.gJFlowRatio    = 8;

cParams.cost         = cost;
cParams.constraint   = constraint;
cParams.initialGuess = x0;
cParams.settings     = s;
cParams.printingPath = true;
problem              = AcademicProblem(cParams);

problem.compute();
xStar = problem.result;

%% Problem 3

close all;
clear;

% Min problem
cost.cF = @(x) -x(1).^2 - x(2).^2 - x(3).^2;
cost.gF = @(x) [-2*x(1); -2*x(2); -2*x(3)];

constraint.cF{1} = @(x) x(1).^2+x(2).^2+2*x(3)-16;
constraint.gF{1} = @(x) [2*x(1); 2*x(2); 2];

constraint.cF{2} = @(x) x(1)+x(2)-4;
constraint.gF{2} = @(x) [1; 1; 0];

% Solution
% x1=2; x2=2; x3=4;


% Setting up
x0               = [0;1;3];
s.type           = "NullSpace";
s.ub             = inf;
s.lb             = 0;
s.maxIter        = 400;
s.constraintCase = {'EQUALITY';'EQUALITY'};
s.etaNorm        = 0.02;
s.gJFlowRatio    = 2;

cParams.cost         = cost;
cParams.constraint   = constraint;
cParams.initialGuess = x0;
cParams.settings     = s;
cParams.printingPath = true;
problem              = AcademicProblem(cParams);

problem.compute();
xStar = problem.result;

%% Problem 4

close all;
clear;

% Min problem
cost.cF = @(x) 2*x(1).^2-3*x(2).^2-2*x(1);
cost.gF = @(x) [4*x(1)-2; -6*x(2)];

constraint.cF{1} = @(x) x(1).^2+x(2).^2-1;
constraint.gF{1} = @(x) [2*x(1); 2*x(2)];

% Solution
% x1=0.2; x2=+-0.9798; l=3;


% Setting up
x0               = [5;2];
s.type           = "NullSpace";
s.ub             = inf;
s.lb             = -inf;
s.maxIter        = 400;
s.constraintCase = {'INEQUALITY'};
s.etaNorm        = 0.01;
s.gJFlowRatio    = 1;

cParams.cost         = cost;
cParams.constraint   = constraint;
cParams.initialGuess = x0;
cParams.settings     = s;
cParams.printingPath = true;
problem              = AcademicProblem(cParams);

problem.compute();
xStar = problem.result;

%% Problem 5

close all;
clear;

% Min problem
cost.cF = @(x) 3*x(1)^2+x(2)^2+x(3)^2-2*x(1)*x(2);
cost.gF = @(x) [6*x(1)-2*x(2); 2*x(2)-2*x(1); 2*x(3)];

constraint.cF{1} = @(x) 2*x(1)^2+x(2)^2+x(3)^2-1;
constraint.gF{1} = @(x) [4*x(1); 2*x(2); 2*x(3)];

% Solution
% x1=x2=x3=0; l=0;


% Setting up
x0               = [1;1;1];
s.type           = "NullSpace";
s.ub             = inf;
s.lb             = -inf;
s.maxIter        = 400;
s.constraintCase = {'INEQUALITY'};
s.etaNorm        = 0.02;
s.gJFlowRatio    = 4;

cParams.cost         = cost;
cParams.constraint   = constraint;
cParams.initialGuess = x0;
cParams.settings     = s;
cParams.printingPath = true;
problem              = AcademicProblem(cParams);

problem.compute();
xStar = problem.result;


%% Test case 3 Paper NullSpaceSLERP

close all;
clear;

% Min problem
cost.cF = @(x) (x(1)+3).^2+x(2).^2;
cost.gF = @(x) [2*(x(1)+3); 2*(x(2))];

constraint.cF{1} = @(x) -x(1).^2+x(2);
constraint.gF{1} = @(x) [-2*x(1); 1];

constraint.cF{2} = @(x) -x(1)-x(2)-2;
constraint.gF{2} = @(x) [-1; -1];

% Solution



% Setting up
x0               = [3;3];
s.type           = "NullSpace";
s.ub             = [4;4];
s.lb             = [-4;1];
s.maxIter        = 400;
s.constraintCase = {'INEQUALITY','INEQUALITY'};
s.etaNorm        = 0.2;
s.gJFlowRatio    = 1;

cParams.cost         = cost;
cParams.constraint   = constraint;
cParams.initialGuess = x0;
cParams.settings     = s;
cParams.printingPath = true;
problem              = AcademicProblem(cParams);

problem.compute();
xStar = problem.result;