% This code generates a dataset of constitutive tensor values for various
% possible combinations of ellipse parameters

close all
clear
clc

% Handle paths
addpath(genpath('Tutorials'))
addpath(genpath('src'))

% Set geometrical parameters
n_variations = 10;
min_semiAxis = 0.01;
max_semiAxis = 0.49;

min_mParam = 1;
max_mParam = 30;

min_nParam = 0.5;
max_nParam = 20;

% Data-file storage
data_filename = 'PauFolder/Datasets/Chomog_superformula_Big.csv';

%% Compute the homogenized tensors



% Compute ellipse sides
gPar.semiHorizontalAxis = 0.45;
gPar.semiVerticalAxis = 0.45;


gPar.m  = 6;
gPar.n1 = 1;
gPar.n2 = 7;
gPar.n3 = 8;

% Find homogenised constitutive tensor
femMicro = SuperformulaDbFEMElasticityMicro(gPar);

% Fetch homogenized constitutive tensor
Chomog_tensor = femMicro.stateProblem.Chomog;

% Plot mesh
femMicro.mesh.plot();

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

adjust_figure_properties(hfig, 12, 30, 0.5);