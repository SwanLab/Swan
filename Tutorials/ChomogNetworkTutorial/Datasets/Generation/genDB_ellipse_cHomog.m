% This code generates a dataset of constitutive tensor values for various
% possible combinations of ellipse parameters

close all
clear
clc

%% Set geometrical parameters
n_variations = 10;
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

%% Data Storage
% Create a table with Chomog_array and Sides_array
%data_table = array2table([Sides_array, Chomog_array], 'VariableNames', {'a', 'b', 'Chomog_00', 'Chomog_01', 'Chomog_02', 'Chomog_10', 'Chomog_11', 'Chomog_12', 'Chomog_20', 'Chomog_21', 'Chomog_22'});
data_table = array2table([Sides_array, Chomog_array(:, [1, 2, 5, 9])], 'VariableNames', {'a', 'b', 'Chomog_00', 'Chomog_01', 'Chomog_11', 'Chomog_22'});

% Save the table as a CSV file
writetable(data_table, data_filename);

%% Data plotting

close all

% Reshape the a and b parameters
var_a = reshape(Sides_array(:, 1), [n_variations, n_variations]);
var_b = reshape(Sides_array(:, 2), [n_variations, n_variations]);

% Plot every Chomog variable (tensor position)
hfig = figure;
tiledlayout(3, 3)
for i = 1:9

    % Reshape the data into a matrix with positions dependant on a and b
    var_C = reshape(Chomog_array(:, i), [n_variations, n_variations]);
    
    % Plot the data
    nexttile
    surf(var_a, var_b, var_C)
    %shading interp
    colormap winter

    % Find the indices of the constitutive tensor (must be a simpler way)
    [C_comp_j, C_comp_i] = ind2sub([3, 3], i);

    % Indicate ellipse parameters - axis
    xlabel('Sx')
    ylabel('Sy')

    % Indicate component of constitutive tensor
    title(['Component ', num2str(C_comp_i), ', ', num2str(C_comp_j),' of constitutive tensor'])
end

adjust_figure_properties(hfig, 12, 30, 0.7);