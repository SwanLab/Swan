clear all;
close all;

BeamOpt = EulerBeamOptimizer();
val = BeamOpt.value;
cost = val(end);

test = cost-1.16
if abs(test) < 0.1
    disp('Test Passed')
else 
    disp('Test Failed')
end