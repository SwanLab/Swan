function buildSidesRotationDensitySurrogate(inFile,outFile,mU)
if nargin < 1 || isempty(inFile)
    inFile = 'latticeSidesRotationNN_Adam_64x3_cholesky.mat';
end
if nargin < 2 || isempty(outFile)
    outFile = 'latticeSidesRotationDensityNN_Adam_64x3_cholesky.mat';
end
if nargin < 3 || isempty(mU)
    mU = 0.99;
end
v = load(inFile);
if ~isequal(v.parameterNames(1:2),{'m1','m2'})
    error('buildSidesRotationDensitySurrogate:NotSides','The first two parameters of %s are not m1 and m2.',inFile);
end
problemNN      = v.problemNN;
surrogate      = DensityAspectSurrogate(v.surrogate,mU);
parameterNames = [{'rho','t'},v.parameterNames(3:end)];
settings       = v.settings;
report         = v.report;
save(outFile,'problemNN','surrogate','parameterNames','settings','report','mU');
rng(1);
n  = 2000;
h  = 1e-6;
X  = [0.03+0.95*rand(n,1),rand(n,1),pi*(rand(n,1)-0.5)];
[~,dY] = surrogate.computeOutputValuesAndGradient(X);
err    = zeros(numel(parameterNames),1);
for i = 1:numel(parameterNames)
    Xp = X; Xp(:,i) = Xp(:,i) + h;
    Xm = X; Xm(:,i) = Xm(:,i) - h;
    fd = (surrogate.computeOutputValues(Xp) - surrogate.computeOutputValues(Xm))/(2*h);
    an = reshape(dY(i,:,:),size(fd,2),[])';
    err(i) = norm(fd-an,'fro')/max(norm(an,'fro'),eps);
end
disp(table(err,'RowNames',parameterNames(:),'VariableNames',{'gradientError'}))
[m1,m2] = surrogate.computeSides([0.3;0.3;0.3],[0;0.5;1]);
fprintf('At rho = 0.3: t = 0 gives (m1,m2) = (%.3f,%.3f), t = 0.5 gives (%.3f,%.3f), t = 1 gives (%.3f,%.3f)\n',m1(1),m2(1),m1(2),m2(2),m1(3),m2(3));
fprintf('Saved %s with parameters %s\n',outFile,strjoin(parameterNames,', '));
end