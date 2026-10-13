clc
clear
close all

rng(1)

rho = linspace(0,1,15);

aMin = 0;
aMax = 0.98;

bMin = 0.3;
bMax = 3;

nGeomTrain = 50;
nGeomTest = 15;

nBoundary = 4;
nGeomLHS = nGeomTrain-nBoundary;

lhsTrain = lhsdesign(nGeomLHS,2);
lhsTest = lhsdesign(nGeomTest,2);

aTrain = ...
    aMin + lhsTrain(:,1)*(aMax-aMin);

bTrain = ...
    bMin + lhsTrain(:,2)*(bMax-bMin);

geometryTrainInterior = ...
    [aTrain,bTrain];

geometryBoundary = [...
    aMin,bMin;
    aMin,bMax;
    aMax,bMin;
    aMax,bMax];

geometryTrain = [...
    geometryTrainInterior;
    geometryBoundary];

aTest = ...
    aMin + lhsTest(:,1)*(aMax-aMin);

bTest = ...
    bMin + lhsTest(:,2)*(bMax-bMin);

geometryTest = ...
    [aTest,bTest];

fprintf('\n')
fprintf('============================================\n')
fprintf('        LATTICE NN DATASET GENERATION\n')
fprintf('============================================\n')

fprintf('\nNumber of rho values: %d\n',length(rho))
fprintf('Training geometries:  %d\n',size(geometryTrain,1))
fprintf('Test geometries:      %d\n',size(geometryTest,1))

fprintf('\nTraining samples: %d\n',...
    length(rho)*size(geometryTrain,1))

fprintf('Test samples:     %d\n',...
    length(rho)*size(geometryTest,1))

fprintf('\nParameter ranges\n')
fprintf('rho: [%g, %g]\n',min(rho),max(rho))
fprintf('a:   [%g, %g]\n',aMin,aMax)
fprintf('b:   [%g, %g]\n',bMin,bMax)

figure

scatter(...
    geometryTrainInterior(:,1),...
    geometryTrainInterior(:,2),...
    40,'filled')

hold on

scatter(...
    geometryBoundary(:,1),...
    geometryBoundary(:,2),...
    70,'filled')

scatter(...
    geometryTest(:,1),...
    geometryTest(:,2),...
    50)

xlabel('a')
ylabel('b')

legend(...
    'Training LHS',...
    'Training boundary',...
    'Test LHS',...
    'Location','best')

title('LHS sampling of lattice geometries')

grid on

xlim([aMin,aMax])
ylim([bMin,bMax])

drawnow

lattice = TutorialLatticeNN(false);

fprintf('\n')
fprintf('Computing training dataset...\n')

tic

trainDataset = ...
    lattice.computeDataset(...
    rho,geometryTrain);

trainTime = toc;

fprintf(...
    'Training dataset completed in %.2f s\n',...
    trainTime)

fprintf('\n')
fprintf('Computing test dataset...\n')

tic

testDataset = ...
    lattice.computeDataset(...
    rho,geometryTest);

testTime = toc;

fprintf(...
    'Test dataset completed in %.2f s\n',...
    testTime)

trainFileName = ...
    'LatticeNNTrainData.csv';

testFileName = ...
    'LatticeNNTestData.csv';

trainGeometryFileName = ...
    'LatticeNNTrainGeometries.csv';

testGeometryFileName = ...
    'LatticeNNTestGeometries.csv';

writematrix(...
    trainDataset,...
    trainFileName)

writematrix(...
    testDataset,...
    testFileName)

writematrix(...
    geometryTrain,...
    trainGeometryFileName)

writematrix(...
    geometryTest,...
    testGeometryFileName)

fprintf('\n')
fprintf('============================================\n')
fprintf('              DATASET SUMMARY\n')
fprintf('============================================\n')

fprintf('\nTraining dataset size\n')
disp(size(trainDataset))

fprintf('Test dataset size\n')
disp(size(testDataset))

fprintf('Training geometry size\n')
disp(size(geometryTrain))

fprintf('Test geometry size\n')
disp(size(geometryTest))

fprintf('\nTraining input ranges\n')

disp([...
    min(trainDataset(:,1:3),[],1);
    max(trainDataset(:,1:3),[],1)])

fprintf('\nTest input ranges\n')

disp([...
    min(testDataset(:,1:3),[],1);
    max(testDataset(:,1:3),[],1)])

fprintf('\nFiles created:\n')
fprintf('%s\n',trainFileName)
fprintf('%s\n',testFileName)
fprintf('%s\n',trainGeometryFileName)
fprintf('%s\n',testGeometryFileName)

fprintf('\nTotal homogenizations: %d\n',...
    size(trainDataset,1)+size(testDataset,1))

fprintf('Total execution time: %.2f s\n',...
    trainTime+testTime)