clear all;
close all;

load('initValue.mat')
s.value = value;
BeamOpt = EulerBeamOptimizer(s);
val = BeamOpt.value;
cost = val(end);

test = cost-51.5
if abs(test) < 1
    disp('Test Passed')
else 
    disp('Test Failed')
end