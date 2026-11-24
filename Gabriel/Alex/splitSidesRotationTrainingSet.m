function split = splitSidesRotationTrainingSet(file,testFraction,seed,prefix)
d      = load(file);
m1     = d.data(:,d.columns == "m1");
m2     = d.data(:,d.columns == "m2");
theta  = d.data(:,d.columns == "theta");
isTwin = theta >= min(theta) + d.period/2 - 1e-8;
key    = [m1,m2,theta];
key(isTwin,:) = [m2(isTwin),m1(isTwin),theta(isTwin)-d.period/2];
[cells,~,cellId] = unique(round(key,8),'rows');
sides      = cells(:,1:2);
isInterior = all(sides > min(sides,[],1) & sides < max(sides,[],1),2);
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
fprintf('%d samples, %d different lattices, %d interior: %d samples for test, %d for training\n',numel(theta),size(cells,1),numel(candidates),split.nTest,split.nTrain);
end

function saveSubset(d,keep,file)
subset            = d;
subset.data       = d.data(keep,:);
subset.rhoNominal = d.rhoNominal(keep);
subset.xi         = d.xi(keep);
subset.tensors    = d.tensors(:,:,:,:,keep);
subset.source     = d.source(keep);
save(file,'-struct','subset');
end