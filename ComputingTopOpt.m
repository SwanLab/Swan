function ComputingTopOpt
close all



s.testName = 'test_cantilever2';%''testJose';
s.testName = 'test_cantilever_nullspace';
s.testName = 'PerimeterAsConstraint';
k.penaltyProv = 40;
k.trustProv   = 100;
t = TopOptComputer(s);
t.compute(k);
end