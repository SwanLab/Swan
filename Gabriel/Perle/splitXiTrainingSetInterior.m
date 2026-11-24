function split = splitXiTrainingSetInterior(file,testFraction,seed,prefix)
if nargin < 1 || isempty(file)
    file = fullfile('latticeXiDataset','dataset.mat');
end
if nargin < 2 || isempty(testFraction)
    testFraction = 0.15;
end
if nargin < 3 || isempty(seed)
    seed = 1;
end
if nargin < 4 || isempty(prefix)
    prefix = 'latticeXi';
end
d          = load(file);
key        = orbitKey(d);
[geometries,~,geometryId] = unique(key,'rows');
% isInterior = abs(geometries(:,1)) < max(abs(geometries(:,1))) & ...
%     geometries(:,2) > min(geometries(:,2)) & geometries(:,2) < max(geometries(:,2));
isInterior = all(geometries > min(geometries,[],1) & geometries < max(geometries,[],1),2);
candidates = find(isInterior);
nTest      = round(testFraction*size(geometries,1));
rng(seed);
testGeometries = candidates(randperm(numel(candidates),nTest));
isTest         = ismember(geometryId,testGeometries);
split.trainFile      = [prefix,'Train.mat'];
split.testFile       = [prefix,'Test.mat'];
split.testGeometries = geometries(testGeometries,:);
split.nTrain         = nnz(~isTest);
split.nTest          = nnz(isTest);
saveSubset(d,~isTest,split.trainFile);
saveSubset(d,isTest,split.testFile);
fprintf('%d geometry orbits (a, t): %d for test, %d for training\n',size(geometries,1),nTest,size(geometries,1)-nTest);
fprintf('Samples: %d train -> %s, %d test -> %s\n',split.nTrain,split.trainFile,split.nTest,split.testFile);
end
function key = orbitKey(d)
transformation = BoostTransformation();
key = round([d.data(:,d.columns == "a"),d.t],10);
for name = d.settings.dataset.symmetries
    symmetry = transformation.getSymmetry(name{1});
    key = min(key,round(symmetry.mapParams(key),10));
end
end
function saveSubset(d,keep,file)
subset.data       = d.data(keep,:);
subset.columns    = d.columns;
subset.tensors    = d.tensors(:,:,:,:,keep);
subset.rhoNominal = d.rhoNominal(keep);
subset.t          = d.t(keep);
subset.source     = d.source(keep);
subset.settings   = d.settings;
save(file,'-struct','subset');
end