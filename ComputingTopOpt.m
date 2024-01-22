clear;

% LS:      1e-3
% Density: 1.25e-3

testNames = ["LSTotIso";"LSRelIso";"DensityTotAni";"DensityRelAni";"LSTotAni";"LSRelAni"];
testaG    = [1e-3;1e-3;1.25e-3;1.25e-3;1e-3;1e-3];
scAngle   = [45;45;85;85;85;85];

for i = 1:length(testNames)
    test = char(testNames(i));
    aG   = testaG(i);
    beta = scAngle(i);
    runTopOpt(test,aG,beta);
end

function runTopOpt(setName,aG,beta)

close all;
s.testName = setName;
s.aG       = aG;
s.beta     = beta;
t = TopOptComputer(s);
t.compute();
p1Params.fValues = t.computation.designVariable.fun.fValues;
p1Params.mesh    = t.computation.designVariable.mesh;
Result           = P1Function(p1Params);
Result.print(['APReview/',setName],'Paraview');
saveas(gcf,['APReview/',setName,'.fig']);

end