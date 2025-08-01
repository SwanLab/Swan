close all
s.meshType = 'Metamaterial';

s.bcProp.nSteps = 76;
s.bcProp.maxVal = 1;
s.bcProp.type   = 'DisplacementTractionX';

s.matProp.mu     = 1;
s.matProp.lambda = 1;

s.monitoring.set       = true;
s.monitoring.printInfo = true;
s.monitoring.printFile = false;
s.monitoring.fileNameOut = 'NeoElastic';

s.tolerance = 1e-12;
s.maxIter   = 100;

tester = TestingHyperelasticity(s);
outputData = tester.compute();
