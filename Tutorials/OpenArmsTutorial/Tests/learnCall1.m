% Test LearnableVariables constructor (calls computeInitialTheta)

% Define neural net architecture
params.neuronsPerLayer = [2, 3, 1];
params.nLayers = numel(params.neuronsPerLayer);

% DO NOT define thetavec — let Julia initialize it!

% Call the Julia constructor
output = callJuliaClass('LearnableVariables', 'LearnableVars', params);

% Display thetavec to check weight values
disp('thetavec:');
disp(output.thetavec);

params.thetavec = output.thetavec;

% Optional: reshape to W and b and display
output2 = callJuliaClass('LearnableVariables', 'reshapeInLayerForm', params);
disp('Weight matrices (W):');
for i = 1:length(output2.W)
    disp(['Layer ', num2str(i), ':']);
    disp(output2.W{i});
end

disp('Bias vectors (b):');
for i = 1:length(output2.b)
    disp(['Layer ', num2str(i), ':']);
    disp(output2.b{i});
end
