clear;

testNames = ["LSPure";"DensityPure"];
testaG    = [0.01;0.01];

for i = 1:length(testNames)
    test = char(testNames(i));
    aG   = testaG(i);
    runTopOpt(test,aG);
end

function runTopOpt(setName,aG)

close all;
s.testName = setName;
s.aG       = aG;
t = TopOptComputer(s);
t.compute();
p1Params.fValues = t.computation.designVariable.fun.fValues;
p1Params.mesh    = t.computation.designVariable.mesh;
Result           = P1Function(p1Params);
Result.print(['APReview/',setName],'Paraview');
saveas(gcf,['APReview/',setName,'.fig']);

end