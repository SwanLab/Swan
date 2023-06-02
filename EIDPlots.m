% Eroded, intermediate and dilated plots

% Inputs:
data.resultsName = 'example';
data.beta = 5;
data.eta  = 0.45;




load([data.resultsName,'.mat']); % .mat file containing "designVariable" from t.computation
printResult('Dilated',designVariable,data);
printResult('Intermediate',designVariable,data);
printResult('Eroded',designVariable,data);

function printResult(type,designVariable,data)

resultsName = data.resultsName;
beta = data.beta;
switch type
    case 'Dilated'
        eta = 0.5 - data.eta;
    case 'Intermediate'
        eta = 0.5;
    case 'Eroded'
        eta = 0.5 + data.eta;
end
rho = designVariable.value(1:end-1);

s.mesh = designVariable.mesh;
s.filterType = 'Filter&Project';
s.designVarType = 'Density';
s.quadratureOrder = 'LINEAR';
s.femSettings = [];
s.femSettings.scale = 'MACRO';
s.femSettings.beta = beta;
s.femSettings.eta  = eta;
Threshold = Filter.create(s);
rhoI = Threshold.getP0fromP1(rho);

p.fValues = mean(rhoI,2);
p.mesh    = s.mesh;
p.filename = [resultsName,type];
p.type = 'GiD';
rhoIFun   = P0Function(p);
rhoIFun.print(p);

end