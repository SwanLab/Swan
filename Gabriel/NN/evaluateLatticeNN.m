function report = evaluateLatticeNN(problemNN,file,s)
if nargin < 3 || isempty(s)
    s = latticeNNSettings();
end
d     = load(file);
X     = d.data(:,s.data.xFeatures);
Y     = d.data(:,s.data.yFeatures);
names = cellstr(d.columns(s.data.yFeatures));
[Ypred,dY] = problemNN.computeOutputValuesAndGradient(X);
report.file       = file;
report.global     = relativeError(Ypred,Y);
report.components = componentErrors(Ypred,Y,names);
report.bins       = arrayfun(@(b) binErrors(Ypred,Y,X,b),s.evaluation.bins,'UniformOutput',false);
report.jacobian   = jacobianSummary(dY,s.parameterNames);
printReport(report);
end

function e = relativeError(Ypred,Y)
e = norm(Ypred-Y,'fro')/max(norm(Y,'fro'),eps);
end

function t = componentErrors(Ypred,Y,names)
err        = Ypred - Y;
relative   = vecnorm(err)./max(vecnorm(Y),eps);
meanAbs    = mean(abs(err));
maxAbs     = max(abs(err));
t = table(relative',meanAbs',maxAbs','RowNames',names,'VariableNames',{'relative','meanAbs','maxAbs'});
end

function t = binErrors(Ypred,Y,X,bin)
values = bin.map(X);
id     = discretize(values,bin.edges);
nBins  = numel(bin.edges) - 1;
lower  = bin.edges(1:end-1)';
upper  = bin.edges(2:end)';
nSamples = accumarray(id(~isnan(id)),1,[nBins 1]);
relative = arrayfun(@(k) relativeError(Ypred(id==k,:),Y(id==k,:)),(1:nBins)');
t = table(lower,upper,nSamples,relative);
t.Properties.Description = bin.name;
end

function j = jacobianSummary(dY,parameterNames)
j.isFinite = all(isfinite(dY(:)));
magnitude  = squeeze(vecnorm(dY,2,2));
j.table    = table(mean(magnitude,2),max(magnitude,[],2),'RowNames',parameterNames,'VariableNames',{'meanNorm','maxNorm'});
end

function printReport(r)
fprintf('\n%s\n',r.file);
fprintf('Global relative error: %.4e\n\n',r.global);
disp(r.components);
for k = 1:numel(r.bins)
    fprintf('Relative error by %s\n',r.bins{k}.Properties.Description);
    disp(r.bins{k});
end
fprintf('Jacobian finite: %d\n',r.jacobian.isFinite);
disp(r.jacobian.table);
end