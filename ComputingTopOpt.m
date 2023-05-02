function ComputingTopOpt

% Note: You can use FEMInputWriter to create benchmarking tests!
% Note: Use gid to create harder tests!

s.testName = 'test_cantilever_IPM';%''testJose';
t = TopOptComputer(s);
t.compute();

end