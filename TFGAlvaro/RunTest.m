s.testName = 'testAlvaro';
t = TopOptComputer(s);
t.compute();



% With the following lines you obtain the result for the last iteration
% (example: design variable with GiD. Test other results and also ParaView!)
p1Params.fValues = t.computation.designVariable.value;
p1Params.mesh    = t.computation.designVariable.mesh;
Result           = P1Function(p1Params);



ResultsName = 'CantileverExample';
type   = 'Paraview'; % GiD/Paraview
Result.print(['TFGAlvaro/Results/',ResultsName],type);