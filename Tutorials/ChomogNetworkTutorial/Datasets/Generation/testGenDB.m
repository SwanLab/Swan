% This code generates a dataset of constitutive tensor values for various
% possible combinations of ellipse parameters

close all
clear
clc

%% Set geometrical parameters
n_variations = 1;
min_semiAxis = 0.01;
max_semiAxis = 0.49;

% Data-file storage
data_filename = 'Tutorials/ChomogNetworkTutorial/Datasets/DB_ellipse_cHomog_testing.csv';

%% Compute the homogenized tensors

% Initialize arrays of ellipse parameters, and results arrays
lengths_array = linspace(min_semiAxis, max_semiAxis, n_variations);
Chomog_array = zeros(length(lengths_array)^2, 9);
Sides_array = zeros(length(lengths_array)^2, 2);

% Loop through different combinations of ellipse sides
counter = 0;
for a = 1:length(lengths_array)
    for b = 1:length(lengths_array)

        % Compute ellipse sides
        gPar.xSide = lengths_array(a);
        gPar.ySide = lengths_array(b);

        % Find homogenised constitutive tensor
        femMicro = EllipseDbFEMElasticityMicro(gPar);

        % Find storage index
        w = b + (a - 1) * length(lengths_array);

        % Fetch homogenized constitutive tensor
        Chomog_mdt = femMicro.stateProblem.Chomog;
        Chomog_tensor = tensorToVoigt2D(Chomog_mdt);

        % Store values of interest
        Chomog_array(w, :) = [Chomog_tensor(1, :), Chomog_tensor(2, :), Chomog_tensor(3, :)];
        Sides_array(w, :) = [gPar.xSide, gPar.ySide];

        counter = counter + 1;
        disp(['Progress: ', num2str(100*counter/length(lengths_array)^2), ' %']);

    end
end