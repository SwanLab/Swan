close all
s.meshType = '';

s.bcProp.nSteps = 30;
s.bcProp.maxVal = -5;
s.bcProp.type   = 'Displacement';

s.matProp.mu     = 1;
s.matProp.lambda = 1;

s.monitoring.set       = true;
s.monitoring.printInfo = true;
s.monitoring.printFile = true;
s.monitoring.fileNameOut = 'NeoElastic';

s.tolerance = 1e-12;
s.maxIter   = 100;

tester = TestingHyperelasticity(s);
outputData = tester.compute();
