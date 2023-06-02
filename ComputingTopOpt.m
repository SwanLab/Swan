function ComputingTopOpt

close all

% Note: You can use FEMInputWriter to create benchmarking tests!
% Note: Use gid to create harder tests!

s.testName = 'test_nullspace';
t = TopOptComputer(s);
t.compute();


% With the following lines you obtain the result for the last iteration
% (example: design variable with GiD. Test other results and also ParaView!)
p1Params.fValues = t.computation.designVariable.value;
p1Params.mesh    = t.computation.designVariable.mesh;
Result           = P1Function(p1Params);
c.type = 'GiD';
c.filename = [s.testName,'_LastIter'];
Result.print(c);

physProblem = t.computation.cost.shapeFunctions{1,1}.getPhysicalProblems();
p.filename = 'testFailureFunInit';
p.type = 'GiD';
physProblem{1,1}.failureFun.print(p);

p.filename = 'testFailureFunStressInit';
physProblem{1,1}.stressFun.print(p);

p.filename = 'testFailureFunStrainInit';
physProblem{1,1}.strainFun.print(p);

p.filename = 'testFailureFunDisplacementInit';
physProblem{1,1}.uFun.print(p);

end