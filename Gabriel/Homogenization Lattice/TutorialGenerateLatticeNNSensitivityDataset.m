clc
clear
close all

rng(2)

rho = linspace(0,1,15);

aMin = 0;
aMax = 0.98;

bMin = 0.3;
bMax = 3;

hRelative = 5e-3;

hA = hRelative*(aMax-aMin);
hB = hRelative*(bMax-bMin);

nGeomSensitivity = 20;

aSampleMin = aMin+hA;
aSampleMax = aMax-hA;

bSampleMin = bMin+hB;
bSampleMax = bMax-hB;

lhsSensitivity = ...
    lhsdesign(nGeomSensitivity,2);

aSensitivity = ...
    aSampleMin + ...
    lhsSensitivity(:,1)*...
    (aSampleMax-aSampleMin);

bSensitivity = ...
    bSampleMin + ...
    lhsSensitivity(:,2)*...
    (bSampleMax-bSampleMin);

geometrySensitivity = ...
    [aSensitivity,bSensitivity];

geometryAPlus = ...
    [aSensitivity+hA,bSensitivity];

geometryAMinus = ...
    [aSensitivity-hA,bSensitivity];

geometryBPlus = ...
    [aSensitivity,bSensitivity+hB];

geometryBMinus = ...
    [aSensitivity,bSensitivity-hB];

fprintf('\n')
fprintf('============================================\n')
fprintf('      LATTICE NN SENSITIVITY DATASET\n')
fprintf('============================================\n')

fprintf('\nNumber of rho values:       %d\n',...
    length(rho))

fprintf('Sensitivity geometries:     %d\n',...
    nGeomSensitivity)

fprintf('Sensitivity samples:        %d\n',...
    length(rho)*nGeomSensitivity)

fprintf('\nhA = %.12e\n',hA)
fprintf('hB = %.12e\n',hB)

fprintf('\nParameter ranges\n')
fprintf('rho: [%g, %g]\n',min(rho),max(rho))
fprintf('a:   [%g, %g]\n',aMin,aMax)
fprintf('b:   [%g, %g]\n',bMin,bMax)

figure

scatter(...
    geometrySensitivity(:,1),...
    geometrySensitivity(:,2),...
    50,'filled')

xlabel('a')
ylabel('b')

title('LHS sampling for sensitivity data')

grid on

xlim([aMin,aMax])
ylim([bMin,bMax])

drawnow

lattice = TutorialLatticeNN(false);

fprintf('\nComputing reference tensor...\n')

tic

baseDataset = ...
    lattice.computeDataset(...
    rho,geometrySensitivity);

baseTime = toc;

fprintf('Reference tensor completed in %.2f s\n',...
    baseTime)

fprintf('\nComputing a+h dataset...\n')

tic

aPlusDataset = ...
    lattice.computeDataset(...
    rho,geometryAPlus);

aPlusTime = toc;

fprintf('a+h completed in %.2f s\n',...
    aPlusTime)

fprintf('\nComputing a-h dataset...\n')

tic

aMinusDataset = ...
    lattice.computeDataset(...
    rho,geometryAMinus);

aMinusTime = toc;

fprintf('a-h completed in %.2f s\n',...
    aMinusTime)

fprintf('\nComputing b+h dataset...\n')

tic

bPlusDataset = ...
    lattice.computeDataset(...
    rho,geometryBPlus);

bPlusTime = toc;

fprintf('b+h completed in %.2f s\n',...
    bPlusTime)

fprintf('\nComputing b-h dataset...\n')

tic

bMinusDataset = ...
    lattice.computeDataset(...
    rho,geometryBMinus);

bMinusTime = toc;

fprintf('b-h completed in %.2f s\n',...
    bMinusTime)

C = baseDataset(:,4:9);

dCda = ...
    (aPlusDataset(:,4:9) - ...
     aMinusDataset(:,4:9))/(2*hA);

dCdb = ...
    (bPlusDataset(:,4:9) - ...
     bMinusDataset(:,4:9))/(2*hB);

sensitivityDataset = [...
    baseDataset(:,1:3),...
    C,...
    dCda,...
    dCdb];

fileName = ...
    'LatticeNNSensitivityData.csv';

geometryFileName = ...
    'LatticeNNSensitivityGeometries.csv';

writematrix(...
    sensitivityDataset,...
    fileName)

writematrix(...
    geometrySensitivity,...
    geometryFileName)

fprintf('\n')
fprintf('============================================\n')
fprintf('          SENSITIVITY DATA SUMMARY\n')
fprintf('============================================\n')

fprintf('\nDataset size\n')
disp(size(sensitivityDataset))

fprintf('Geometry size\n')
disp(size(geometrySensitivity))

fprintf('\nMaximum absolute C\n')
disp(max(abs(C),[],1))

fprintf('\nMaximum absolute dC/da\n')
disp(max(abs(dCda),[],1))

fprintf('\nMaximum absolute dC/db\n')
disp(max(abs(dCdb),[],1))

fprintf('\nMean norm C      = %.12e\n',...
    mean(vecnorm(C,2,2)))

fprintf('Mean norm dC/da = %.12e\n',...
    mean(vecnorm(dCda,2,2)))

fprintf('Mean norm dC/db = %.12e\n',...
    mean(vecnorm(dCdb,2,2)))

totalTime = ...
    baseTime + ...
    aPlusTime + ...
    aMinusTime + ...
    bPlusTime + ...
    bMinusTime;

fprintf('\nFiles created:\n')
fprintf('%s\n',fileName)
fprintf('%s\n',geometryFileName)

fprintf('\nTotal FEM evaluations: %d\n',...
    5*size(baseDataset,1))

fprintf('Total execution time: %.2f s\n',...
    totalTime)