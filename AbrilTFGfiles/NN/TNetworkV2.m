% This script is the same one as the 1st version but it creates a loop and
% trains all the columns at once saving it all compactly in one variable to
% upload on the preconditioner.

clc;
clear;
close all;

%% Initialization of hyperparameters
pol_deg         = 1;
testratio       = 30;
lambda          = 0.0;
learningRate    = 0.1;
hiddenLayers    = [24 40 60 60 40 24 12];
 

%% INITIALIZATION 
% Store dataset file name
s.fileName = 'DataT.csv';

% Load model parameters
s.polynomialOrder = pol_deg;
s.testRatio       = testratio;
s.networkParams.hiddenLayers    = hiddenLayers;
s.optimizerParams.learningRate  = learningRate;
s.optimizerParams.maxEpochs = 1000; % 1000 is the best option, but we use 10 to pass the tutorial quickly
s.costParams.lambda             = lambda;
s.costParams.costType           = 'L2';

s.networkParams.HUtype = 'ReLU';
s.networkParams.OUtype = 'linear';


% Select the model's features
s.xFeatures = 1:3;


%% Initialization of variables to save
T_NN=cell(1,8);

%% Loop for the 8 coarse modes

for j=8:8
    switch j
        case 1
            s.yFeatures = [4,5];     
            TitleName="Cost Function T1";
        case 2         
            s.yFeatures = [6,7];      
            TitleName="Cost Function T2";
        case 3         
            s.yFeatures = [8,9];     
            TitleName="Cost Function T3";
        case 4   
            s.yFeatures = [10,11];    
            TitleName="Cost Function T4";
        case 5  
            s.yFeatures = [12,13];    
            TitleName="Cost Function T5";
        case 6   
            s.yFeatures = [14,15];   
            TitleName="Cost Function T6";
        case 7
            s.yFeatures = [16,17];    
            TitleName="Cost Function T7";
        case 8  
            s.yFeatures = [18,19];    
            TitleName="Cost Function T8";
    end
    
    
    % Load data
    data   = Data(s);
    s.data = data;
    
    % Train the model
    T_NN{j} = OptimizationProblemNN(s);
    T_NN{j}.solve();
    T_NN{j}.plotCostFnc();
    title(TitleName);

end

FileName=fullfile('AbrilTFGfiles','NN',"T_NN2.mat");
save(FileName, "T_NN");

