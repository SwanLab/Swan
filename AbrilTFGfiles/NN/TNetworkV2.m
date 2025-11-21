% This script is the same one as the 1st version but it creates a loop and
% trains all the columns at once saving it all compactly in one variable to
% upload on the preconditioner.

clc;
clear;
close all;

%% Initialization of hyperparameters
pol_deg         = 1;
testratio       = 20;
lambda          = 0.0;
learningRate    = 0.01;
hiddenLayers    = [4,16,16,16];
 

%% INITIALIZATION 
% Store dataset file name
s.fileName = 'DataT.csv';

% Load model parameters
s.polynomialOrder = pol_deg;
s.testRatio       = testratio;
s.networkParams.hiddenLayers    = hiddenLayers;
s.optimizerParams.learningRate  = learningRate;
s.optimizerParams.maxEpochs = 2000; % 1000 is the best option, but we use 10 to pass the tutorial quickly
s.costParams.lambda             = lambda;
s.costParams.costType           = 'L2';

s.networkParams.HUtype = 'ReLU';
s.networkParams.OUtype = 'linear';


% Select the model's features
s.xFeatures = 1:3;


%% Initialization of variables to save
MSETrain=zeros(1,8);
T_NN=cell(1,8);
comparison=cell(1,8);

%% Loop for the 8 coarse modes

for j=1:8
    switch j
        case 1
            s.yFeatures = [4,5];      %T1
            TitleName="Cost Function T1";
            FileName="T1.mat";
        case 2         
            s.yFeatures = [6,7];      %T2
            TitleName="Cost Function T2";
            FileName="T2.mat";
        case 3         
            s.yFeatures = [8,9];      %T3
            TitleName="Cost Function T3";
            FileName="T3.mat";
        case 4   
            s.yFeatures = [10,11];    %T4
            TitleName="Cost Function T4";
            FileName="T4.mat";
        case 5  
            s.yFeatures = [12,13];    %T5
            TitleName="Cost Function T5";
            FileName="T5.mat";
        case 6   
            s.yFeatures = [14,15];    %T6
            TitleName="Cost Function T6";
            FileName="T6.mat";
        case 7
            s.yFeatures = [16,17];    %T7
            TitleName="Cost Function T7";
            FileName="T7.mat";
        case 8  
            s.yFeatures = [18,19];    %T8
            TitleName="Cost Function T8";
            FileName="T8.mat";
    end
    
    
    % Load data
    data   = Data(s);
    s.data = data;
    
    % Train the model
    T_NN{1,j} = OptimizationProblemNN(s);
    T_NN{1,j}.solve();
    
   % Mirar si funciona 
   % FileName=fullfile('AbrilTFGfiles','NN',FileName);
   % save(FileName, "T_NN{1,j}");

    T_NN{1,j}.plotCostFnc();
    title(TitleName);
    
    MSETrain(1,j) = immse(T_NN{1,j}.computeOutputValues(data.Xtrain), data.Ytrain);
    
    
    % --- COMPARISON ---
    
    % Load dataset from specified path
    filePath = fullfile('AbrilTFGfiles', s.fileName);
    tempData = readmatrix(filePath);
    
    comparison{1,j}.real = tempData(:, s.yFeatures);
    comparison{1,j}.predicted = zeros(size(comparison{1,j}.real));
    
    for i = 1:size(comparison{1,j}.real,1)
        comparison{1,j}.predicted(i, :) = T_NN{1,j}.computeOutputValues(tempData(i, s.xFeatures));
    end
    
    comparison{1,j}.difference = comparison{1,j}.real-comparison{1,j}.predicted;

end




FileName=fullfile('AbrilTFGfiles','NN',"T_NN.mat");
    save(FileName, "T_NN","MSETrain","comparison");

