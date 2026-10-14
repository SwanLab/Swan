clc
clear
close all

rng(1)

s.trainFileName = ...
    'LatticeNNTrainData.csv';

s.testFileName = ...
    'LatticeNNTestData.csv';

s.polynomialOrder = 3;

s.xFeatures = 1:3;
s.yFeatures = 4:9;

data = Data(s);

networkParams.hiddenLayers = ...
    [40 72 150 72];

networkParams.HUtype = ...
    'SoftPlus';

networkParams.OUtype = ...
    'linear';

networkParams.data = data;

network = Network(networkParams);

designVariable = ...
    network.getLearnableVariables();

lossParams.network = network;

lossParams.designVariable = ...
    designVariable;

lossParams.data = data;

lossParams.costType = 'L2';

loss = LossFunctional(lossParams);

theta0 = designVariable.thetavec;

[J,g] = ...
    loss.computeFunctionAndGradient(theta0);

fprintf('\n')
fprintf('============================================\n')
fprintf('      SOFTPLUS BACKPROP VALIDATION\n')
fprintf('============================================\n')

fprintf('\nInitial cost = %.12e\n',J)

fprintf('Number of parameters = %d\n',...
    length(theta0))

fprintf('Size theta\n')
disp(size(theta0))

fprintf('Size gradient\n')
disp(size(g))

fprintf('Gradient norm = %.12e\n',...
    norm(g))

fprintf('Maximum gradient = %.12e\n',...
    max(abs(g)))

fprintf('\n')
fprintf('============================================\n')
fprintf('       DIRECTIONAL DERIVATIVE TEST\n')
fprintf('============================================\n')

nDirections = 5;

hValues = [...
    1e-4,...
    1e-5,...
    1e-6];

for iDirection = 1:nDirections

    direction = ...
        randn(size(theta0));

    direction = ...
        direction/norm(direction);

    dAnalytical = ...
        sum(g.*direction);

    fprintf('\nDirection %d\n',iDirection)

    fprintf(...
        'Analytical derivative = %.12e\n',...
        dAnalytical)

    for iH = 1:length(hValues)

        h = hValues(iH);

        thetaP = ...
            theta0 + h*direction;

        thetaM = ...
            theta0 - h*direction;

        [Jp,~] = ...
            loss.computeFunctionAndGradient(...
            thetaP);

        [Jm,~] = ...
            loss.computeFunctionAndGradient(...
            thetaM);

        dFD = ...
            (Jp-Jm)/(2*h);

        absError = ...
            abs(dAnalytical-dFD);

        relError = ...
            absError/...
            max(abs(dFD),eps);

        fprintf(...
            ['h = %.1e   ',...
             'FD = %.12e   ',...
             'Abs = %.12e   ',...
             'Rel = %.12e\n'],...
            h,...
            dFD,...
            absError,...
            relError)

    end

end

fprintf('\n')
fprintf('============================================\n')
fprintf('       PARAMETER-BY-PARAMETER TEST\n')
fprintf('============================================\n')

nPL = network.neuronsPerLayer;

testIndices = [];

parameterNames = {};

last = 1;

for iLayer = 2:length(nPL)

    nW = ...
        nPL(iLayer-1)*nPL(iLayer);

    nB = ...
        nPL(iLayer);

    wStart = last;

    wEnd = ...
        wStart+nW-1;

    bStart = ...
        wEnd+1;

    bEnd = ...
        bStart+nB-1;

    wMiddle = ...
        floor((wStart+wEnd)/2);

    indicesLayer = [...
        wStart,...
        wMiddle,...
        wEnd,...
        bStart,...
        bEnd];

    namesLayer = {...
        sprintf('W%d first',iLayer-1),...
        sprintf('W%d middle',iLayer-1),...
        sprintf('W%d last',iLayer-1),...
        sprintf('b%d first',iLayer-1),...
        sprintf('b%d last',iLayer-1)};

    testIndices = ...
        [testIndices,indicesLayer];

    parameterNames = ...
        [parameterNames,namesLayer];

    last = bEnd+1;

end

h = 1e-6;

fprintf('\n')
fprintf(['Parameter         Index       ',...
         'Analytical          FD', ...
         '                 Abs error', ...
         '          Rel error\n'])

for iTest = 1:length(testIndices)

    index = testIndices(iTest);

    thetaP = theta0;
    thetaM = theta0;

    thetaP(index) = ...
        thetaP(index)+h;

    thetaM(index) = ...
        thetaM(index)-h;

    [Jp,~] = ...
        loss.computeFunctionAndGradient(...
        thetaP);

    [Jm,~] = ...
        loss.computeFunctionAndGradient(...
        thetaM);

    gFD = ...
        (Jp-Jm)/(2*h);

    gAnalytical = ...
        g(index);

    absError = ...
        abs(gAnalytical-gFD);

    relError = ...
        absError/...
        max(abs(gFD),eps);

    fprintf(...
        '%-15s %8d   % .12e   % .12e   %.12e   %.12e\n',...
        parameterNames{iTest},...
        index,...
        gAnalytical,...
        gFD,...
        absError,...
        relError)

end

designVariable.thetavec = theta0;

[Jfinal,gFinal] = ...
    loss.computeFunctionAndGradient(theta0);

fprintf('\n')
fprintf('============================================\n')
fprintf('               FINAL CHECK\n')
fprintf('============================================\n')

fprintf('\nInitial cost = %.12e\n',J)
fprintf('Final cost   = %.12e\n',Jfinal)

fprintf(...
    'Theta restoration error = %.12e\n',...
    norm(designVariable.thetavec-theta0))

fprintf(...
    'Gradient restoration error = %.12e\n',...
    norm(gFinal-g))

fprintf('\nBackprop validation completed.\n')