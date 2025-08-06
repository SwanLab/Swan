% This code generates a dataset of constitutive tensor values for various
% possible combinations of ellipse parameters

close all
clear
clc
% Set geometrical parameters
n_variations = 100;
min_semiAxis = 0.0;
max_semiAxis = 0.5;
xSlope = 0.2;
ySlope = 3.9;
nIntercept = -0.3;

% Data-file storage
data_filename = 'Tutorials/ChomogNetworkTutorial/Datasets/DB_plane_slope_testin.csv';

%% Compute the homogenized tensors

% Initialize arrays of ellipse parameters, and results arrays
lengths_array = linspace(min_semiAxis, max_semiAxis, n_variations);
Chomog_array = zeros(length(lengths_array)^2, 1);
Sides_array = zeros(length(lengths_array)^2, 2);

% Loop through different combinations of ellipse sides
counter = 0;
for a = 1:length(lengths_array)
    for b = 1:length(lengths_array)

        % Find inputs
        xSide = lengths_array(a);
        ySide = lengths_array(b);

        % Find storage index
        w = b + (a - 1) * length(lengths_array);

        % Fetch homogenized constitutive tensor
        Yout = nIntercept + xSlope * xSide + ySlope * ySide;

        % Store values of interest
        Chomog_array(w, :) = [Yout];
        Sides_array(w, :) = [xSide, ySide];

        counter = counter + 1;
        disp(['Progress: ', num2str(100*counter/length(lengths_array)^2), ' %']);

    end
end

%% Data Storage
% Create a table with Chomog_array and Sides_array
data_table = array2table([Sides_array, Chomog_array(:)], 'VariableNames', {'a', 'b', 'yOut'});

% Save the table as a CSV file
writetable(data_table, data_filename);

%% Data plotting

close all

% Reshape the a and b parameters
var_a = reshape(Sides_array(:, 1), [n_variations, n_variations]);
var_b = reshape(Sides_array(:, 2), [n_variations, n_variations]);

% Plot every Chomog variable (tensor position)
hfig = figure;
% Reshape the data into a matrix with positions dependant on a and b
var_C = reshape(Chomog_array(:), [n_variations, n_variations]);

% Plot the data
surf(var_a, var_b, var_C)
%shading interp
colormap winter

% Find the indices of the constitutive tensor (must be a simpler way)

% Indicate ellipse parameters - axis
xlabel('Sx')
ylabel('Sy')