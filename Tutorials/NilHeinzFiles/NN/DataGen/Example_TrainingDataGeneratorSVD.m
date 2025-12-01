% EXAMPLE_TRAININGDATAGENERATORSVD Ejemplo de uso de TrainingDataGenerator con SVD
% 
% Este ejemplo muestra cómo usar TrainingDataGenerator para generar datos
% de entrenamiento con el enfoque SVD optimizado.

close all;
clear;
clc;

%% Configuración
radii = 0.5*(0.1:0.05:0.9);  % Radios para el conjunto de entrenamiento
fprintf('=== Ejemplo: Generación de Datos con SVD ===\n\n');

%% Crear generador de datos con parámetros SVD personalizados
% Parámetros opcionales:
%   'svdTolerance': Tolerancia para truncamiento (default: 1e-6)
%   'svdEnergyRatio': Ratio de energía para truncamiento (default: 0.99)
dataGen = TrainingDataGenerator(radii, ...
    'svdTolerance', 1e-6, ...
    'svdEnergyRatio', 0.99);

%% Generar datos (incluyendo SVD)
% computeSVD = true para calcular SVD, false para omitir
fprintf('Generando datos de entrenamiento...\n');
dataGen.generateData(true);  % true = calcular SVD

%% Exportar datos SVD
fprintf('\nExportando datos SVD...\n');
dataGen.exportSVDToCSV(outputDir);  % Exporta VTrainingData.csv
dataGen.exportSVDToMAT(outputDir, 'SVD_Results.mat');  % Exporta U, S, V

%% Obtener información sobre SVD
k = dataGen.getNumberOfModes();
fprintf('\n=== Resumen SVD ===\n');
fprintf('Número de modos retenidos: %d\n', k);
fprintf('Reducción de dimensionalidad: De %d parámetros → %d parámetros\n', ...
    8 * size(dataGen.allResults{1}.u, 1), k);

%% Verificar error de reconstrucción para algunos radios
fprintf('\n=== Verificación de Reconstrucción ===\n');
testIndices = [1, round(length(radii)/2), length(radii)];
for idx = testIndices
    err = dataGen.computeReconstructionError(idx);
    fprintf('Radio r=%.3f (índice %d): Error relativo = %.2e\n', ...
        radii(idx), idx, err);
end

%% Visualización opcional: Valores singulares
[U, S, V, k] = dataGen.getSVDResults();
sigma = diag(S);

figure('Position', [100, 100, 1200, 400]);

% Subplot 1: Valores singulares (escala logarítmica)
subplot(1, 3, 1);
semilogy(1:length(sigma), sigma, 'o-', 'LineWidth', 1.5, 'MarkerSize', 6);
hold on;
semilogy([1, k], [sigma(k), sigma(k)], 'r--', 'LineWidth', 1.5);
xlabel('Modo i');
ylabel('Valor singular σ_i');
title('Valores Singulares');
legend('σ_i', sprintf('Modo %d (truncado)', k), 'Location', 'best');
grid on;

% Subplot 2: Energía acumulada
subplot(1, 3, 2);
energy = cumsum(sigma.^2) / sum(sigma.^2);
plot(1:length(energy), energy * 100, 'o-', 'LineWidth', 1.5, 'MarkerSize', 6);
hold on;
plot([1, k], [energy(k)*100, energy(k)*100], 'r--', 'LineWidth', 1.5);
xlabel('Número de modos');
ylabel('Energía acumulada (%)');
title('Energía Acumulada');
legend('Energía', sprintf('Modo %d (%.1f%%)', k, energy(k)*100), 'Location', 'best');
grid on;
ylim([0, 100]);

% Subplot 3: Primeros modos paramétricos V
subplot(1, 3, 3);
nModesToPlot = min(5, k);
colors = lines(nModesToPlot);
hold on;
for i = 1:nModesToPlot
    plot(radii, V(:, i), '-', 'LineWidth', 1.5, 'Color', colors(i,:), ...
        'DisplayName', sprintf('V_%d', i));
end
xlabel('Radio r');
ylabel('Coeficiente V_i');
title('Modos Paramétricos (primeros 5)');
legend('Location', 'best');
grid on;

sgtitle('Análisis SVD de Matrices T', 'FontSize', 14, 'FontWeight', 'bold');

fprintf('\n=== Proceso completado ===\n');
fprintf('Archivos generados en: %s\n', outputDir);
fprintf('  - TTrainingData.csv\n');
fprintf('  - KcoarseTrainingData.csv\n');
fprintf('  - McoarseTrainingData.csv\n');
fprintf('  - VTrainingData.csv (datos SVD para entrenamiento)\n');
fprintf('  - SVD_Results.mat (U, S, V para reconstrucción)\n');

