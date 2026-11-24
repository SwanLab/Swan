function split = splitParallelogramTrainingSet(file,testFraction,seed,prefix)
d      = load(file);
names  = ["m1","m2","a","s","r"];
params = zeros(size(d.data,1),numel(names));
for i = 1:numel(names)
    params(:,i) = d.data(:,d.columns == names(i));
end
own    = round(params(:,[4,5,3,1,2]),6);
twin   = TurnBoostStretchTransformation().computeTwin(params);
twin   = round(twin(:,[4,5,3,1,2]),6);
delta  = own - twin;
[~,j]  = max(delta ~= 0,[],2);
first  = delta(sub2ind(size(delta),(1:size(delta,1))',j));
key    = own;
key(first > 0,:) = twin(first > 0,:);
[cells,~,cellId] = unique(key,'rows');
sides      = cells(:,4:5);
isInterior = all(sides > min(sides,[],1) & sides < max(sides,[],1),2);
candidates = find(isInterior);
nTest      = round(testFraction*size(cells,1));
rng(seed);
testCells  = candidates(randperm(numel(candidates),nTest));
isTest     = ismember(cellId,testCells);
split.trainFile = [prefix,'Train.mat'];
split.testFile  = [prefix,'Test.mat'];
split.nTrain    = nnz(~isTest);
split.nTest     = nnz(isTest);
saveSubset(d,~isTest,split.trainFile);
saveSubset(d,isTest,split.testFile);
fprintf('%d samples, %d different lattices, %d interior: %d samples for test, %d for training\n',size(params,1),size(cells,1),numel(candidates),split.nTest,split.nTrain);
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