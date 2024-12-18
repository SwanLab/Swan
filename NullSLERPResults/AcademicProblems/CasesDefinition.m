%% Case 1
close all;
clear;

% Min problem
cost.cF = @(x) 0.3*x(1)+x(2);
cost.gF = @(x) [0.3; 1];

constraint.cF{1} = @(x) 1./x(1)-x(2);
constraint.gF{1} = @(x) [-1./((x(1)).^2); -1];

constraint.cF{2} = @(x) x(1)+x(2)-3;
constraint.gF{2} = @(x) [1; 1];

% Setting up
x0               = [1.5;2.25];
s.type           = "NullSpace";
s.ub             = [1.5;2];
s.lb             = [0.8;1];
s.maxIter        = 400;
s.constraintCase = {'INEQUALITY','INEQUALITY'};
s.etaNorm        = 1;
s.tauMax         = 1;
s.gJFlowRatio    = 1; % Only this

cParams.cost         = cost;
cParams.constraint   = constraint;
cParams.initialGuess = x0;
cParams.settings     = s;
cParams.printingPath = true;
problem              = AcademicProblem(cParams);

problem.compute();
xStar = problem.result;



%% Case 2
close all;
clear;

% Min problem
cost.cF = @(x) (x(1)-2).^2+(x(2)-2).^2;
cost.gF = @(x) [2*(x(1)-2); 2*(x(2)-2)];

constraint.cF{1} = @(x) 1./x(1)-x(2);
constraint.gF{1} = @(x) [-1./((x(1)).^2); -1];

constraint.cF{2} = @(x) x(1)+x(2)-3;
constraint.gF{2} = @(x) [1; 1];

% Setting up
x0               = [1.5;2.25];
s.type           = "NullSpace";
s.ub             = [2;1.6];
s.lb             = [1.6;1];
s.maxIter        = 15;
s.constraintCase = {'INEQUALITY','INEQUALITY'};
s.etaNorm        = 1;
s.tauMax         = 0.5;
s.gJFlowRatio    = 5;% 0.2, 1, 5. Explicar necesidad de etaMin para 0.2 y que sólo aplicará aquí

cParams.cost         = cost;
cParams.constraint   = constraint;
cParams.initialGuess = x0;
cParams.settings     = s;
cParams.printingPath = true;
problem              = AcademicProblem(cParams);

problem.compute();
xStar = problem.result;



%% Case 3
close all;
clear;

% Min problem
cost.cF = @(x) (x(1)+3).^2+x(2).^2;
cost.gF = @(x) [2*(x(1)+3); 2*(x(2))];

constraint.cF{1} = @(x) -x(1).^2+x(2);
constraint.gF{1} = @(x) [-2*x(1); 1];

constraint.cF{2} = @(x) -x(1)-x(2)-2;
constraint.gF{2} = @(x) [-1; -1];

% Setting up
x0               = [3;3];
s.type           = "NullSpace";
s.ub             = [4;4];
s.lb             = [-4;1];
s.maxIter        = 400;
s.constraintCase = {'INEQUALITY','INEQUALITY'};
s.etaNorm        = 1;
s.tauMax         = 0.5;
s.gJFlowRatio    = 1; % Only this

cParams.cost         = cost;
cParams.constraint   = constraint;
cParams.initialGuess = x0;
cParams.settings     = s;
cParams.printingPath = true;
problem              = AcademicProblem(cParams);

problem.compute();
xStar = problem.result;