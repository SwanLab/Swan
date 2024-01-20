clear;

testNames = ["LSPure";"DensityPure"];
testaJ    = [0.3;0.6];

for i = 1:length(testNames)
    test = char(testNames(i));
    aJ   = testaJ(i);
    runTopOpt(test,aJ);
end

function runTopOpt(setName,aJ)

close all;
s.testName = setName;
s.aJ       = aJ;
t = TopOptComputer(s);
t.compute();
p1Params.fValues = t.computation.designVariable.fun.fValues;
p1Params.mesh    = t.computation.designVariable.mesh;
Result           = P1Function(p1Params);
Result.print(['APReview/',setName],'Paraview');
saveas(gcf,['APReview/',setName,'.fig']);

end