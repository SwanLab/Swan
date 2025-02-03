%% Consumption Optimization

clc;
clear;
close all;

%% Initialization of hyperparameters
pol_deg         = 1;                
testratio       = 20;              
lambda          = 0;                
learningRate    = 0.1;              
hiddenLayers    = [3,2,1,2,3];      

%% INITIALIZATION 
s.fileName = '/OpenArms/resultadov3.csv';          
s.polynomialOrder = pol_deg;
s.testRatio       = testratio;
s.xFeatures = [1, 2, 3, 4, 5, 6, 7];                
s.yFeatures = [8];                 
data = Data(s);

s.networkParams.costType = 'L2';
s.networkParams.HUtype = 'ReLU';
s.networkParams.OUtype = 'linear';
s.networkParams.hiddenLayers    = hiddenLayers;
s.data                          = data;
s.optimizerParams.learningRate  = learningRate;
s.costParams.lambda             = lambda;

opt = OptimizationProblem(s);
opt.solve();

%% TEST DATA
[Xtest, Ytest] = opt.getTestData();
network = opt.getNetwork();
Ypred = zeros(size(Ytest));

for i = 1:size(Xtest, 1)
    Ypred(i) = network.forwardprop(Xtest(i, :), Ytest(i, :));
end

%% TABLE
difference = Ytest - Ypred;
input_data = array2table(Xtest);
input_data.Properties.VariableNames = {'Rpm', 'Windy(deg)', 'Windy(m/s)', 'Speed (m/s)', 'Yaw', 'Pitch', 'Roll'};

output_data = table(Ytest, Ypred, difference);
output_data.Properties.VariableNames = {'Cons. real', 'Cons. prediction', 'Difference'};

result_table = [input_data, output_data];

disp('Results table');
disp(result_table);

%% MSE calculation
mse = mean((Ypred - Ytest).^2);  
disp(['MSE: ', num2str(mse)]);


