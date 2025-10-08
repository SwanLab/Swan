function genDB_stiffness_matrix()
% Genera un dataset de entrenamiento para redes neuronales que predicen
% matrices de rigidez 8x8 basadas en un parámetro de entrada
%
% OUTPUT: Guarda un archivo CSV con:
%         - Columna 1: parámetro de entrada (ej: densidad, espesor, etc.)
%         - Columnas 2-65: elementos de la matriz de rigidez 8x8 aplanada

clc;
clear;
close all;

%% Parámetros de generación
n_samples = 1000;  % Número de muestras
input_range = [0.1, 2.0];  % Rango del parámetro de entrada
output_file = 'DB_stiffness_matrix.csv';

%% Generar parámetros de entrada
input_params = linspace(input_range(1), input_range(2), n_samples)';

%% Pre-allocar matriz de salida
stiffness_data = zeros(n_samples, 65);  % 1 entrada + 64 salidas
stiffness_data(:, 1) = input_params;

%% Generar matrices de rigidez para cada muestra
fprintf('Generando %d muestras de matrices de rigidez 8x8...\n', n_samples);

for i = 1:n_samples
    % Generar matriz de rigidez basada en el parámetro de entrada
    K = generateStiffnessMatrix(input_params(i));
    
    % Aplanar la matriz a vector
    K_vector = StiffnessMatrixHandler.matrixToVector(K);
    
    % Almacenar en la matriz de datos
    stiffness_data(i, 2:end) = K_vector;
    
    % Mostrar progreso
    if mod(i, 100) == 0
        fprintf('Progreso: %d/%d (%.1f%%)\n', i, n_samples, 100*i/n_samples);
    end
end

%% Crear encabezados para el CSV
headers = {'input_param'};
for i = 1:8
    for j = 1:8
        headers{end+1} = sprintf('K_%d_%d', i, j);
    end
end

%% Guardar datos
output_path = fullfile('Tutorials', 'ChomogNetworkTutorial', 'Datasets', output_file);
writematrix([headers; num2cell(stiffness_data)], output_path);

fprintf('Dataset guardado en: %s\n', output_path);
fprintf('Dimensiones: %d muestras x %d características\n', size(stiffness_data, 1), size(stiffness_data, 2));

%% Visualizar algunas muestras
visualizeSamples(stiffness_data, input_params);

end

function K = generateStiffnessMatrix(input_param)
% Genera una matriz de rigidez 8x8 basada en un parámetro de entrada
% INPUT: input_param - parámetro de entrada (ej: densidad, espesor)
% OUTPUT: K - matriz de rigidez 8x8 válida

% Crear una matriz base que depende del parámetro de entrada
% Esta es una función de ejemplo - deberías reemplazarla con tu modelo físico

% Matriz base con estructura típica de elementos finitos
K_base = [
    2, -1,  0,  0, -1,  0,  0,  0;
   -1,  2, -1,  0,  0, -1,  0,  0;
    0, -1,  2, -1,  0,  0, -1,  0;
    0,  0, -1,  2,  0,  0,  0, -1;
   -1,  0,  0,  0,  2, -1,  0,  0;
    0, -1,  0,  0, -1,  2, -1,  0;
    0,  0, -1,  0,  0, -1,  2, -1;
    0,  0,  0, -1,  0,  0, -1,  2
];

% Escalar la matriz base por el parámetro de entrada
K_scaled = input_param * K_base;

% Agregar variaciones aleatorias para simular diferentes configuraciones
noise_scale = 0.1 * input_param;
K_noisy = K_scaled + noise_scale * randn(8, 8);

% Asegurar que sea una matriz de rigidez válida
K = StiffnessMatrixHandler.makeValidStiffnessMatrix(K_noisy);

end

function visualizeSamples(stiffness_data, input_params)
% Visualiza algunas muestras del dataset generado

figure;
tiledlayout(2, 3);

% Seleccionar algunas muestras para visualizar
sample_indices = [1, 250, 500, 750, 1000];

for i = 1:min(5, length(sample_indices))
    idx = sample_indices(i);
    
    nexttile;
    
    % Extraer matriz de rigidez
    K_vector = stiffness_data(idx, 2:end);
    K = StiffnessMatrixHandler.vectorToMatrix(K_vector);
    
    % Visualizar matriz
    imagesc(K);
    colorbar;
    title(sprintf('Input: %.2f', input_params(idx)));
    xlabel('DOF');
    ylabel('DOF');
end

sgtitle('Ejemplos de Matrices de Rigidez Generadas');

% Plot de algunas propiedades vs parámetro de entrada
figure;
tiledlayout(2, 2);

% Diagonal principal
nexttile;
plot(input_params, stiffness_data(:, 2), 'b-', 'LineWidth', 2);
xlabel('Parámetro de Entrada');
ylabel('K(1,1)');
title('Elemento Diagonal K(1,1)');

% Elemento fuera de la diagonal
nexttile;
plot(input_params, stiffness_data(:, 3), 'r-', 'LineWidth', 2);
xlabel('Parámetro de Entrada');
ylabel('K(1,2)');
title('Elemento K(1,2)');

% Traza de la matriz
nexttile;
trace_values = zeros(size(input_params));
for i = 1:length(input_params)
    K_vector = stiffness_data(i, 2:end);
    K = StiffnessMatrixHandler.vectorToMatrix(K_vector);
    trace_values(i) = trace(K);
end
plot(input_params, trace_values, 'g-', 'LineWidth', 2);
xlabel('Parámetro de Entrada');
ylabel('Traza de K');
title('Traza de la Matriz de Rigidez');

% Número de condición
nexttile;
condition_numbers = zeros(size(input_params));
for i = 1:length(input_params)
    K_vector = stiffness_data(i, 2:end);
    K = StiffnessMatrixHandler.vectorToMatrix(K_vector);
    [condition_numbers(i), ~, ~] = StiffnessMatrixHandler.analyzeMatrix(K);
end
semilogy(input_params, condition_numbers, 'm-', 'LineWidth', 2);
xlabel('Parámetro de Entrada');
ylabel('Número de Condición');
title('Número de Condición de K');

sgtitle('Propiedades de las Matrices de Rigidez vs Parámetro de Entrada');

end
