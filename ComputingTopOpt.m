function ComputingTopOpt

% Note: You can use FEMInputWriter to create benchmarking tests!
% Note: Use gid to create harder tests!

s.testName = 'test_gripping';
t = TopOptComputer(s);
t.compute();

end