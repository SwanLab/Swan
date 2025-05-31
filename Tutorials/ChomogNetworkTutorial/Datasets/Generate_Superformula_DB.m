% This code generates a dataset of constitutive tensor values for various
% possible combinations of ellipse parameters

close all
clear
clc

% Handle paths
addpath(genpath('Tutorials'))
addpath(genpath('src'))

% Set geometrical parameters

max_radius = 0.245;
min_radius = 0.1;

min_semiAxis = 0.2;
max_semiAxis = 0.4;
nVar_semiAxis = 1;

min_mParam = 6;
max_mParam = 10;
nVar_mParam = 1 + (max_mParam - min_mParam) / 2;

min_nParam = 2;
max_nParam = 12;
nVar_nParam = 6;

% Data-file storage
data_filename = 'PauFolder/Datasets/Chomog_superformula_Big.csv';

% Compute the homogenized tensors

semiAxisArray = linspace(min_semiAxis, max_semiAxis, nVar_semiAxis);
mParamArray = min_mParam:2:max_mParam;
nParamArray = linspace(min_nParam, max_nParam, nVar_nParam);

%nVar_total = nVar_semiAxis^2 * nVar_mParam * nVar_nParam^3;

nVar_total = 0;
for n_sh = 1:nVar_semiAxis
    for n_mp = 1:nVar_mParam
        for n_n1 = 1:nVar_nParam
            for n_n2 = 1:nVar_nParam
                for n_n3 = 1:nVar_nParam
                    % Fetch superformula parameters
                    gPar.semiVerticalAxis = semiAxisArray(n_sh); % b
                    gPar.m  = mParamArray(n_mp);
                    gPar.n1 = nParamArray(n_n1);
                    gPar.n2 = nParamArray(n_n2);
                    gPar.n3 = nParamArray(n_n3);

                    % S'ajusta a per un b donat per aconseguir un radi desitjat!!
                    gPar.semiHorizontalAxis = gPar.semiVerticalAxis^(gPar.n3/gPar.n2);
                    
                    if fun_super_eval(gPar, min_radius, max_radius)
                        nVar_total = nVar_total + 1;
                    end
                end
            end
        end
    end
end
fun_logRunTime(nVar_total);

Chomog_array = zeros(nVar_total, 4);
Params_array = zeros(nVar_total, 6);

n_counter = 0;

for n_sh = 1:nVar_semiAxis
    for n_mp = 1:nVar_mParam
        for n_n1 = 1:nVar_nParam
            for n_n2 = 1:nVar_nParam
                for n_n3 = 1:nVar_nParam
                    % Fetch superformula parameters
                    gPar.semiVerticalAxis = semiAxisArray(n_sh); % b
                    gPar.m  = mParamArray(n_mp);
                    gPar.n1 = nParamArray(n_n1);
                    gPar.n2 = nParamArray(n_n2);
                    gPar.n3 = nParamArray(n_n3);

                    % S'ajusta a per un b donat per aconseguir un radi desitjat!!
                    gPar.semiHorizontalAxis = gPar.semiVerticalAxis^(gPar.n3/gPar.n2);
                    %gPar.semiHorizontalAxis = find_a_super(gPar, wanted_maxRadius);

                    if fun_super_eval(gPar, min_radius, max_radius)
                        % Update counter
                        n_counter = n_counter + 1;

                        % Find homogenised constitutive tensor
                        femMicro = SuperformulaDbFEMElasticityMicro(gPar);
                        
                        % Fetch homogenized constitutive tensor
                        Chomog_tensor = femMicro.stateProblem.Chomog;

                        % Store data
                        Chomog_array(n_counter, :) = [Chomog_tensor(1, 1), ...
                                                      Chomog_tensor(1, 2), ...
                                                      Chomog_tensor(2, 2), ...
                                                      Chomog_tensor(3, 3)];
                        Params_array(n_counter, :) = cell2mat(struct2cell(gPar))';

                        % Plot mesh
                        figure();
                        femMicro.mesh.plot();

                        % Log generation progress
                        fun_logProgress(n_counter, nVar_total);
                    end

                end
            end
        end
    end
end

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

function is_valid = fun_super_eval(gPar, r_min, r_max)

    is_valid = true;
    tol = 1e-14;

    a = gPar.semiHorizontalAxis;
    b = gPar.semiVerticalAxis;   
    
    m = gPar.m;
    n1 = gPar.n1;
    n2 = gPar.n2;
    n3 = gPar.n3;
    
    periodicity_residue = abs(abs(cos(pi*m/2))^n2+ a ^n2/(b^n3)*abs(sin(pi*m/2))^n3 - 1);
    large_size_residue = max(fun_superMaxRad(a, b, m, n1, n2, n3) - r_max, 0);
    small_size_residue = max(r_min - fun_superMinRad(a, b, m, n1, n2, n3), 0);

    if periodicity_residue > tol || large_size_residue > 0 || small_size_residue > 0
        is_valid = false;
    end

end

function a = find_a_super(gPar, r_target)

    b = gPar.semiVerticalAxis;   

    m = gPar.m;
    n1 = gPar.n1;
    n2 = gPar.n2;
    n3 = gPar.n3;

    obj_fun = @(a) fun_superMaxRad(a, b, m, n1, n2, n3) - r_target;

    try
        % Root-finding using fzero with initial guess interval
        a = fzero(obj_fun, [1e-16, 10]);

        % Post-processing: discard solution if it's too small or inconsistent
        if a < 1e-6 || fun_superMinRad(a, b, m, n1, n2, n3) < 1e-2
            a = -1;
        end

    catch
        % If fzero fails (e.g., no sign change or non-convergence), return failure indicator
        a = -1;
    end

end


function rad = fun_superform(phi, a, b, m, n1, n2, n3)

    rad = (abs(cos(m.*phi./4)./a).^n2 + abs(sin(m.*phi./4)./b).^n3).^(-1/n1);

end

function rad = fun_superMaxRad(a, b, m, n1, n2, n3)

    rad = max(fun_superform(linspace(0, 2*pi, 1000), a, b, m, n1, n2, n3));

end

function rad = fun_superMinRad(a, b, m, n1, n2, n3)

    rad = min(fun_superform(linspace(0, 2*pi, 1000), a, b, m, n1, n2, n3));

end

function plot_superform(gPar)

    a = gPar.semiHorizontalAxis;
    b = gPar.semiVerticalAxis;   
    
    m = gPar.m;
    n1 = gPar.n1;
    n2 = gPar.n2;
    n3 = gPar.n3;

    phi_vec = linspace(0, 2*pi, 1000);
    rad_vec = fun_superform(phi_vec, a, b, m, n1, n2, n3);
    x_vec = rad_vec .* cos(phi_vec);
    y_vec = rad_vec .* sin(phi_vec);
    figure()
    plot(x_vec, y_vec)
    axis equal

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
    
    currentStatus = round(i / nPoints * 100, 10);
    previousStatus = round((i - 1) / nPoints * 100, 10);
    
    if i == 1
        fprintf('Simulation started.\n');
    elseif i == nPoints
        fprintf('Simulation ended.\n');
    else
        crossed = (previousStatus < keyPoints) & (currentStatus >= keyPoints);
        if any(crossed)
            fprintf('Simulation progress: %.1f %%\n', keyPoints(find(crossed, 1)));
        end
    end

end
