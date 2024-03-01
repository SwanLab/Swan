function ComputingTopOpt

% Note: You can use FEMInputWriter to create benchmarking tests!
% Note: Use gid to create harder tests!

s.testName = 'test_nullspace';
t = TopOptComputer(s);
t.compute();


% With the following lines you obtain the result for the last iteration
% (example: design variable with GiD. Test other results and also ParaView!)
p1Params.fValues = t.computation.designVariable.value;
p1Params.mesh    = t.computation.designVariable.mesh;
p1Params.order   = 'P1';
Result           = LagrangianFunction(p1Params);
c.type = 'GiD';
c.filename = [s.testName,'_LastIter'];
Result.print(c);

end