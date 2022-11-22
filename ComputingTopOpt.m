function ComputingTopOpt

% Conclusion: f and line search in another .m file

clear
clc
close all

k=[];

s.testName = 'PerimeterAsConstraint';%''testJose';
t = TopOptComputer(s);
t.compute(k);

end