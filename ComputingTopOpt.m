function ComputingTopOpt

close all
clear

% Note: You can use FEMInputWriter to create benchmarking tests!
% Note: Use gid to create harder tests!

s.testName = 'test_arturo';
t = TopOptComputer(s);
t.compute();

end