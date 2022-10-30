clear all;
close all;

BeamOpt = EulerBeamOptimizer();
cost = BeamOpt.cost;

test = cost-52
if abs(test) < 1
    disp('Test Passed')
else 
    disp('Test Failed')
end