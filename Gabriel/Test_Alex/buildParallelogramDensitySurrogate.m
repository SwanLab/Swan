function buildParallelogramDensitySurrogate(inFile,outFile,mU)
if nargin < 1 || isempty(inFile)
    inFile = 'latticeParallelogramNN_Adam_128x3_cholesky.mat';
end
if nargin < 2 || isempty(outFile)
    outFile = strrep(inFile,'ParallelogramNN','ParallelogramDensityNN');
end
if nargin < 3 || isempty(mU)
    mU = 0.99;
end
v = load(inFile);
if ~isequal(v.parameterNames(1:2),{'m1','m2'})
    error('buildParallelogramDensitySurrogate:NotSides','The first two parameters of %s are not m1 and m2.',inFile);
end
problemNN      = v.problemNN;
surrogate      = DensityAspectSurrogate(v.surrogate,mU);
parameterNames = [{'rho','t'},v.parameterNames(3:end)];
settings       = v.settings;
report         = v.report;
save(outFile,'problemNN','surrogate','parameterNames','settings','report','mU');
rng(1);
n     = 2000;
h     = 1e-6;
angle = -40 + 170*rand(n,1);
X     = [0.03+0.95*rand(n,1),rand(n,1),sind(angle)./(sind(angle)+cosd(angle)),0.6*(2*rand(n,1)-1),1.5.^(2*rand(n,1)-1)];
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
fprintf('Saved %s with parameters %s\n',outFile,strjoin(parameterNames,', '));
end