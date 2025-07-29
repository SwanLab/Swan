s.meshType = 'Bending';

s.bcProp.nSteps = 20;
s.bcProp.maxVal = 1;
s.bcProp.type   = 'ForceTractionX';

s.matProp.mu     = 1;
s.matProp.lambda = 1;

s.monitoring.set   = true;
s.monitoring.print = true;

s.tolerance = 1e-6;
s.maxIter   = 100;

tester = TestingHyperelasticity(s);
outputData = tester.compute();
