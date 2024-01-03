s.testName = 'testAlvaro';
t = TopOptComputer(s);
t.compute();

% With the following lines you obtain the result for the last iteration
% (example: design variable with GiD. Test other results and also ParaView!)
Result = t.computation.designVariable.fun;
u      = t.computation.physicalProblem.uFun;

ResultsName = 'Cubechair_newcode';
type   = 'Paraview'; % GiD/Paraview
Result.print(['TFGAlvaro/Results/',ResultsName],type);
u.print(['TFGAlvaro/Results/',ResultsName,'Displacements'],type);