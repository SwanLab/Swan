function [problemNN,report,surrogate] = trainSidesNN(s)
if nargin < 1 || isempty(s)
    s = latticeSidesNNSettings();
end
rng(s.seed);
data = Data(s.data);
p.data            = data;
p.networkParams   = s.networkParams;
p.optimizerParams = s.optimizerParams;
p.costParams      = s.costParams;
problemNN = OptimizationProblemNN(p);
problemNN.solve();
surrogate    = buildSidesSurrogate(problemNN,s);
report.train = evaluateLatticeNN(surrogate,s.reference.trainFileName,s);
report.test  = evaluateLatticeNN(surrogate,s.reference.testFileName,s);
parameterNames = s.parameterNames;
settings       = s;
save(s.outputFile,'problemNN','surrogate','parameterNames','settings','report');
fprintf('Relative error: train %.4e, test %.4e\n',report.train.global,report.test.global);
fprintf('Saved %s\n',s.outputFile);
problemNN.plotCostFnc();
end