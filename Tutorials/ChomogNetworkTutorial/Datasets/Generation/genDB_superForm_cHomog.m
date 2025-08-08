% This code generates a dataset of constitutive tensor values for various
% possible combinations of ellipse parameters

close all
clear
clc

%% Set geometrical parameters

sf = superformula_functionality;

max_radius = 0.4;
min_radius = 0.1;

min_semiAxis = 0.2;
max_semiAxis = 0.4;
nVar_semiAxis = 50;

min_mParam = 8;
max_mParam = 8;
nVar_mParam = 1 + (max_mParam - min_mParam);

min_nParam = 2;
max_nParam = 18;
nVar_nParam = 100;

% Data-file storage
data_filename = 'Tutorials/ChomogNetworkTutorial/Datasets/DB_superForm_cHomog.csv';

%% Compute the homogenized tensors

semiAxisArray = linspace(min_semiAxis, max_semiAxis, nVar_semiAxis);
mParamArray = min_mParam:1:max_mParam;
nParamArray = linspace(min_nParam, max_nParam, nVar_nParam);

%nVar_total = nVar_semiAxis^2 * nVar_mParam * nVar_nParam^3;

nVar_total = 0;
for n_sh = 1:nVar_semiAxis
    for n_mp = 1:nVar_mParam
        for n_n1 = 1:nVar_nParam
            %for n_n2 = 1:nVar_nParam
                %for n_n3 = 1:nVar_nParam
                    % Fetch superformula parameters
                    gPar.semiVerticalAxis = semiAxisArray(n_sh); % b
                    gPar.m  = mParamArray(n_mp);
                    gPar.n1 = nParamArray(n_n1);
                    gPar.n2 = nParamArray(n_n1);
                    gPar.n3 = nParamArray(n_n1);
                    %gPar.n2 = nParamArray(n_n2);
                    %gPar.n3 = nParamArray(n_n3);

                    % S'ajusta a per un b donat per aconseguir un radi desitjat!!
                    gPar.semiHorizontalAxis = gPar.semiVerticalAxis^(gPar.n3/gPar.n2);
                    
                    if sf.evaluate(gPar, min_radius, max_radius)
                        nVar_total = nVar_total + 1;
                    end
                %end
            %end
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
            %for n_n2 = 1:nVar_nParam
                %for n_n3 = 1:nVar_nParam
                    % Fetch superformula parameters
                    gPar.semiVerticalAxis = semiAxisArray(n_sh); % b
                    gPar.m  = mParamArray(n_mp);
                    gPar.n1 = nParamArray(n_n1);
                    gPar.n2 = nParamArray(n_n1);
                    gPar.n3 = nParamArray(n_n1);
                    %gPar.n2 = nParamArray(n_n2);
                    %gPar.n3 = nParamArray(n_n3);

                    % S'ajusta a per un b donat per aconseguir un radi desitjat!!
                    gPar.semiHorizontalAxis = gPar.semiVerticalAxis^(gPar.n3/gPar.n2);
                    %gPar.semiHorizontalAxis = find_a_super(gPar, wanted_maxRadius);

                    if sf.evaluate(gPar, min_radius, max_radius)
                        % Update counter
                        n_counter = n_counter + 1;

                        % Find homogenised constitutive tensor
                        femMicro = SuperformulaDbFEMElasticityMicro(gPar);
                        
                        % Fetch homogenized constitutive tensor
                        Chomog_mdt = femMicro.stateProblem.Chomog;
                        Chomog_tensor = tensor_to_voigt_2D(Chomog_mdt);

                        % Store data
                        Chomog_array(n_counter, :) = [Chomog_tensor(1, 1), ...
                                                      Chomog_tensor(1, 2), ...
                                                      Chomog_tensor(2, 2), ...
                                                      Chomog_tensor(3, 3)];

                        Params_array(n_counter, :) = [gPar.semiHorizontalAxis, ...
                                                      gPar.semiVerticalAxis, ...
                                                      gPar.m, ...
                                                      gPar.n1, ...
                                                      gPar.n2, ...
                                                      gPar.n3];

                        % Plot mesh
                        % close all;
                        % figure();
                        % femMicro.mesh.plot();
                        % plot_superform(gPar);

                        % Log generation progress
                        logProgress(n_counter, nVar_total);
                    end

                %end
            %end
        end
    end
end

% Data Storage
% Create a table with Chomog_array and Sides_array
data_table = array2table([Params_array, Chomog_array], 'VariableNames', {'a', 'b', 'm', 'n1', 'n2', 'n3', 'Chomog_00', 'Chomog_01', 'Chomog_11', 'Chomog_22'});

% Save the table as a CSV file
writetable(data_table, data_filename);


%% Useful functions

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
    
    plot(x_vec, y_vec)
    axis equal

end


function fun_logRunTime(nVar_total)

    runTimeSecs = nVar_total * 3.4; % 3.4 are the secs approx for the superformula FEM problem
    runTimeMins = runTimeSecs / 60;
    runTimeHours = runTimeMins / 60;
    runTimeDays = runTimeHours / 24;
    
    if runTimeMins > 1 && runTimeMins < 60
        fprintf('Expected execution time: %.1f minutes', runTimeMins);
    elseif runTimeHours > 1 && runTimeHours < 24
        fprintf('Expected execution time: %.1f hours', runTimeHours);
    elseif runTimeDays > 1
        fprintf('Expected execution time: %.1f days', runTimeDays);
    else
        fprintf('Expected execution time: %.1f seconds', runTimeSecs);
    end
    fprintf('\n');

end