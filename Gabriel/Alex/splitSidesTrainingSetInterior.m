function split = splitSidesTrainingSetInterior(file,testFraction,seed,prefix)
d          = load(file);
params     = round(d.data(:,ismember(d.columns,["m1","m2","a"])),10);
[cells,~,cellId] = unique(params,'rows');
isInterior = all(cells > min(cells,[],1) & cells < max(cells,[],1),2);
candidates = find(isInterior);
nTest      = round(testFraction*size(cells,1));
rng(seed);
testCells  = candidates(randperm(numel(candidates),nTest));
isTest     = ismember(cellId,testCells);
split.trainFile = [prefix,'Train.mat'];
split.testFile  = [prefix,'Test.mat'];
split.testCells = cells(testCells,:);
split.nTrain    = nnz(~isTest);
split.nTest     = nnz(isTest);
saveSubset(d,~isTest,split.trainFile);
saveSubset(d,isTest,split.testFile);
fprintf('%d cells (m1, m2, a), %d interior: %d for test, %d for training\n',size(cells,1),numel(candidates),split.nTest,split.nTrain);
end

function saveSubset(d,keep,file)
subset.data       = d.data(keep,:);
subset.columns    = d.columns;
subset.rhoNominal = d.rhoNominal(keep);
subset.xi         = d.xi(keep);
subset.tensors    = d.tensors(:,:,:,:,keep);
subset.source     = d.source(keep);
subset.settings   = d.settings;
save(file,'-struct','subset');
end