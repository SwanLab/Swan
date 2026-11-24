function buildHoleSidesSurrogate(inFile,outFile)
if nargin < 1 || isempty(inFile)
    inFile = 'latticeXiNN_Adam_64x3_cholesky.mat';
end
if nargin < 2 || isempty(outFile)
    outFile = 'latticeSidesNN_Adam_64x3_cholesky.mat';
end
v = load(inFile);
if ~isequal(v.parameterNames(:)',{'rho','a','xi'})
    error('buildHoleSidesSurrogate:WrongParameters', ...
        'Expected parameters rho, a, xi in %s and found %s.', ...
        inFile,strjoin(v.parameterNames,', '));
end
problemNN      = v.problemNN;
surrogate      = HoleSidesSurrogate(v.surrogate);
parameterNames = {'m1','m2','a'};
settings       = v.settings;
report         = v.report;
save(outFile,'problemNN','surrogate','parameterNames','settings','report');
fprintf('Saved %s with parameters %s\n',outFile,strjoin(parameterNames,', '));
end