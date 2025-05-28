% This code generates a dataset of constitutive tensor values for various
% possible combinations of ellipse parameters

close all
clear
clc

% Handle paths
addpath(genpath('Tutorials'))
addpath(genpath('src'))

% Set geometrical parameters

% min_semiAxis = 0.3;
% max_semiAxis = 0.499;
% nVar_semiAxis = 1;
% 
% min_mParam = 1;
% max_mParam = 30;
% nVar_mParam = 1;
% 
% min_nParam = 0.5;
% max_nParam = 20;
% nVar_nParam = 2;

min_semiAxis = 0.3;
max_semiAxis = 0.4;
nVar_semiAxis = 1;

min_mParam = 2;
max_mParam = 12;
nVar_mParam = 1 + range([min_mParam, max_mParam]);

min_nParam = 2;
max_nParam = 12;
nVar_nParam = 1;

% Data-file storage
data_filename = 'PauFolder/Datasets/Chomog_superformula_Big.csv';

%% Compute the homogenized tensors

semiAxisArray = linspace(min_semiAxis, max_semiAxis, nVar_semiAxis);
mParamArray = min_mParam:max_mParam;
nParamArray = linspace(min_nParam, max_nParam, nVar_nParam);

nVar_total = nVar_semiAxis^2 * nVar_mParam * nVar_nParam^3;
fun_logRunTime(nVar_total);

Chomog_array = zeros(nVar_total, 4);
Params_array = zeros(nVar_total, 6);

n_counter = 0;

for n_sh = 1:nVar_semiAxis
    for n_rad = 1:nVar_semiAxis
        for n_mp = 1:nVar_mParam
            for n_n1 = 1:nVar_nParam
                for n_n2 = 1:nVar_nParam
                    for n_n3 = 1:nVar_nParam
                        % Fetch superformula parameters
                        gPar.semiHorizontalAxis = semiAxisArray(n_sh);
                        wanted_maxRadius = semiAxisArray(n_rad);   
                        
                        gPar.m  = mParamArray(n_mp);
                        gPar.n1 = nParamArray(n_n1);
                        gPar.n2 = nParamArray(n_n2);
                        gPar.n3 = nParamArray(n_n3);

                        %stockRad = max(fun_superMaxRad(gPar));
                        %% COM AJUSTO B Per un radi donat??????
                        gPar.semiVerticalAxis = wanted_maxRadius^(gPar.n1/gPar.n3);  

                        % Find homogenised constitutive tensor
                        femMicro = SuperformulaDbFEMElasticityMicro(gPar);
                        
                        % Fetch homogenized constitutive tensor
                        Chomog_tensor = femMicro.stateProblem.Chomog;

                        % Store data
                        n_counter = n_counter + 1;
                        Chomog_array(n_counter, :) = [Chomog_tensor(1, 1), ...
                                                      Chomog_tensor(1, 2), ...
                                                      Chomog_tensor(2, 2), ...
                                                      Chomog_tensor(3, 3)];
                        Params_array(n_counter, :) = cell2mat(struct2cell(gPar))';

                        % Log generation progress
                        fun_logProgress(n_counter, nVar_total);

                    end
                end
            end
        end
    end
end

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


%% Random function

function rad = fun_superMaxRad(gPar)

a = gPar.semiHorizontalAxis;
b = gPar.semiVerticalAxis;   

m = gPar.m;
n1 = gPar.n1;
n2 = gPar.n2;
n3 = gPar.n3;

fun_rad = @(phi) (abs(cos(m.*phi./4)./a).^n2 + abs(sin(m.*phi./4)./b).^n3).^(-1/n1);
rad = max(fun_rad(linspace(0, 2*pi, 500)));

end

function fun_logRunTime(nVar_total)

    runTimeSecs = nVar_total * 3.4;
    runTimeMins = runTimeSecs / 60;
    runTimeHours = runTimeMins / 60;
    runTimeDays = runTimeHours / 24;
    
    if runTimeMins > 1 && runTimeMins < 60
        fprintf('Expected execution time: %.1f minutes', runTimeMins);
    elseif runTimeHours > 1 && runTimeHours < 24
        fprintf('Expected execution time: %.1f hours', runTimeMins);
    elseif runTimeDays > 1
        fprintf('Expected execution time: %.1f days', runTimeDays);
    else
        fprintf('Expected execution time: %.1f seconds', runTimeSecs);
    end
    fprintf('\n');

end

function fun_logProgress(i, nPoints)

    increment = 2.5; % Percentage increment
    keyPoints = increment:increment:100;
    
    currentStatus = i / nPoints * 100;
    previouStatus = (i - 1) / nPoints * 100;
    
    if i == 1
        fprintf('Simulation started.\n')
    elseif i == nPoints
        fprintf('Simulation ended.\n')
    else
        currentState = (currentStatus>keyPoints) .* (previouStatus<keyPoints);
        [~, currentIdx] = max(currentState);
        if any(currentState)
            fprintf('Simulation progress: %.f %% \n', keyPoints(currentIdx));
        end
    end

end