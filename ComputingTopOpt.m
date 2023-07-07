function ComputingTopOpt

close all

% Note: You can use FEMInputWriter to create benchmarking tests!
% Note: Use gid to create harder tests!

s.testName = 'test_cantilever_IPM';%''testJose';
t = TopOptComputer(s);
t.compute();

p.mesh = t.computation.designVariable.mesh;
p.fValues = t.computation.designVariable.value;
Result = P1Function(p);
q.filename = 'results';
q.type     = 'GiD';
Result.print(q);

end