function buildXiTSurrogate(inFile,outFile,mU)
if nargin < 1 || isempty(inFile)
    inFile = 'latticeXiNN_Adam_64x3_cholesky.mat';
end
if nargin < 2 || isempty(outFile)
    outFile = 'latticeXiTNN_Adam_64x3_cholesky.mat';
end
if nargin < 3 || isempty(mU)
    mU = latticeXiSettings().hole.mU;
end
v = load(inFile);
names = v.parameterNames;
if ~strcmp(names{end},'xi')
    error('buildXiTSurrogate:NotXi','The last parameter of %s is %s, not xi.',inFile,names{end});
end
problemNN      = v.problemNN;
surrogate      = XiBoundsSurrogate(v.surrogate,mU);
parameterNames = [names(1:end-1),{'t'}];
settings       = v.settings;
settings.hole.mU = mU;
report         = v.report;
save(outFile,'problemNN','surrogate','parameterNames','settings','report');
fprintf('Saved %s with parameters %s\n',outFile,strjoin(parameterNames,', '));
end