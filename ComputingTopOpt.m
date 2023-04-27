function ComputingTopOpt

clear
clc
close all

s.testName = 'test_gripping';
t = TopOptComputer(s);
t.compute();

end