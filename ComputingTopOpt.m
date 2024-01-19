
runTopOpt('LSPure');

function runTopOpt(setName)

close all;
s.testName = setName;
t = TopOptComputer(s);
t.compute();
p1Params.fValues = t.computation.designVariable.fun.fValues;
p1Params.mesh    = t.computation.designVariable.mesh;
Result           = P1Function(p1Params);
Result.print(['APReview/',setName],'Paraview');
saveas(gcf,['APReview/',setName,'.fig']);

end