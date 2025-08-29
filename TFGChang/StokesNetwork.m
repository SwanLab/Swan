clc;
clear;
close all;

%% Preprocess and Save Dataset

OriginalData = readmatrix('E_AoA5_mpt.txt');
data         = Data.preProcessDataset(OriginalData);

save("EAoA5Datas.mat","data");

%% Initialize and Train the NN

for i = 1:1:1

    % Store dataset file name
    s.fileName          = "EDatasNS.mat";
    s.testRatio         = 20;
    s.k                 = i;
    s.polynomialOrder   = 1;

    % Select the model's features
    s.xFeatures = [1, 2, 3, 4];
    s.yFeatures = 7;

    % Load data
    data   = Data(s);
    s.data = data;

    % Initialization

    % Load model parameters
    s.networkParams.hiddenLayers    = s.polynomialOrder * size(s.xFeatures,2)* 9 * ones(1,6);
    s.optimizerParams.learningRate  = 0.05;
    s.costParams.lambda             = 0;
    s.costParams.costType           = 'L2';

    s.networkParams.HUtype = 'ReLU';
    s.networkParams.OUtype = 'linear';

    % Train the model
    opt = OptimizationProblem(s);
    opt.solve();
    opt.plotCostFnc();

    % Compute Error

    MSETrain(i)    = immse(opt.computeOutputValues(data.Xtrain), data.Ytrain);
    MSETest(i)     = immse(opt.computeOutputValues(data.Xtest), data.Ytest);

    disp("MSE:")
    disp(MSETrain)
    disp(MSETest)

end

% Compute Mean Error

disp("------")
disp(mean(MSETrain));
disp(mean(MSETest));

%% Init for Plot

dataset  = readmatrix("E_AoA5_mpt.txt");

normalized = true;
if (normalized == true)
    load("StokesNN.mat");
    EData = readmatrix("EData.txt");
    maxValue        = max(EData(:,end));
    minValue        = min(EData(:,end));
    dataset(:,end)  = (dataset(:,end) - minValue) / (maxValue - minValue);
else
    load("StokesNNO.mat")
end

%% Plot Symmetric Airfoil E Data Comparison

yDataSym = opt.computeOutputValues(dataset(1:16, 1:4));

% Extract Original Data
tSym = dataset(1:16,3);
eSym = dataset(1:16,end);

figure;
plot(tSym, eSym);
hold on;
plot(tSym, yDataSym);
legend('Location', 'best');
title("Aerodynamic efficiency of symmetric airfoil vs thickness (t) - AoA 5º");
xlabel('Thickness (t)');
ylabel('Aerodynamic Efficiency');
legend('Ground Truth', 'Predicted');

figure;
ReError = (yDataSym - eSym) ./  eSym;
plot(tSym, ReError);
title("Relative Error between Ground Truth and Predicted Aerodynamic Efficiency of Symmetric Airfoil vs Thickness (t)");
xlabel('Thickness (t)');
ylabel('Relative Error');

%% Plot Asymetric Airfoil E Surface

% Extract Original Data
mT = dataset(17:end,1);
pT = dataset(17:end,2);
tT = dataset(17:end,3);
eT = dataset(17:end,end);

m = unique(mT);
p = unique(pT);

[M,P] = meshgrid(m,p);

% Thickness Values to Plot
thicknesses = [0.1, 0.2, 0.3, 0.4];

figure;

for i = 1:length(thicknesses)
    ti = thicknesses(i);

    idx = abs(tT - ti) < 1e-3;

    if (i == 1)
        mi = mT(idx);
        pi = pT(idx);
    end

    ei = eT(idx);

    F = scatteredInterpolant(mi, pi, ei, 'linear', 'none');
    E(:,:,i) = F(M, P);

    subplot(2, 2, i);
    surface(M,P,E(:,:,i));
    shading interp;
    title(sprintf('t = %.2f', ti));
    xlabel('Max Camber (m)');
    ylabel('Max Camber Position (p)');
    axis tight;
    colorbar;
end

sgtitle('Ground truth aerodynamic efficiency (E) vs. m and p for different thickness (t) - AoA 5º','FontWeight', 'bold');

% Compute Predicted Outputs y

yData = opt.computeOutputValues(dataset(17:end, 1:4));
figure;

for i = 1:length(thicknesses)
    ti = thicknesses(i);

    idx = abs(tT - ti) < 1e-3;
    ei = yData(idx);

    F = scatteredInterpolant(mi, pi, ei, 'linear', 'none');
    EP(:,:,i) = F(M, P);

    subplot(2, 2, i);
    surface(M,P,EP(:,:,i));
    shading interp;
    title(sprintf('t = %.2f', ti));
    xlabel('Max Camber (m)');
    ylabel('Max Camber Position (p)');
    axis tight;
    colorbar;
end

sgtitle('Predicted aerodynamic efficiency vs. m and p for different thickness (t) - AoA 5º','FontWeight', 'bold');

ReError = abs(EP - E) ./ E;

figure;

for i = 1:length(thicknesses)
    ti = thicknesses(i);

    subplot(2, 2, i);
    surface(M,P,ReError(:,:,i));
    shading interp;
    title(sprintf('t = %.2f', ti));
    xlabel('Max Camber (m)');
    ylabel('Max Camber Position (p)');
    axis tight;
    colorbar;
end

sgtitle('Relative error between prediction and ground truth vs. m and p for different thickness (t) - AoA 5º','FontWeight', 'bold');

%% Plot Merged E Surface

% Compute Predicted Outputs
yData  = opt.computeOutputValues(dataset(17:end, 1:4));

% Extract Original Data
mT = dataset(17:end,1);
pT = dataset(17:end,2);
tT = dataset(17:end,3);
eT = dataset(17:end,end);

m = unique(mT);
p = unique(pT);

[M,P] = meshgrid(m,p);

% Thickness Values to Plot
thicknesses = [0.1, 0.2, 0.3, 0.4];

figure;

for i = 1:length(thicknesses)
    ti = thicknesses(i);

    idx = abs(tT - ti) < 1e-3;

    if (i == 1)
        mi = mT(idx);
        pi = pT(idx);
    end

    eiT = eT(idx);
    eiP = yData(idx);

    F = scatteredInterpolant(mi, pi, eiT, 'linear', 'none');
    E(:,:,i) = F(M, P);

    F = scatteredInterpolant(mi, pi, eiP, 'linear', 'none');
    EP(:,:,i) = F(M, P);

    subplot(2, 2, i);
    surface(M,P,E(:,:,i),'FaceColor', 'blue', 'EdgeColor', 'none', 'DisplayName', 'Ground Truth');
    surface(M,P,EP(:,:,i),'FaceColor', [1, 0.5, 0], 'EdgeColor', 'none', 'DisplayName', 'Predicted');
    alpha(0.7);
    view(125,20)

    legend('Location', 'best');
    %shading interp;
    title(sprintf('t = %.2f', ti));
    xlabel('Max Camber (m)');
    ylabel('Max Camber Position (p)');
    axis tight;
    %colorbar;
end

sgtitle('Ground truth and predicted efficiency vs. m and p for different thickness (t)','FontWeight', 'bold');
