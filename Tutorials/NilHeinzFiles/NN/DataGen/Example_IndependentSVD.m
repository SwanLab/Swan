% EXAMPLE_INDEPENDENTSVD Ejemplo de uso independiente de SVD
% 
% Este ejemplo demuestra cómo generar datos sin SVD primero, y luego
% calcular SVD de forma independiente.

close all;
clear;
clc;

fprintf('=== Ejemplo: Generación Independiente de Datos y SVD ===\n\n');

%% Configuración
radii = 0.1:0.05:0.9;  % Radios para el conjunto de entrenamiento

%% ================================================================
%% FLUJO 1: Generar datos SIN SVD primero
%% ================================================================
fprintf('--- FLUJO 1: Generación sin SVD ---\n');

% Crear generador de datos
dataGen = TrainingDataGenerator(radii);

% Generar datos SIN SVD
fprintf('\nGenerando datos de entrenamiento (sin SVD)...\n');
dataGen.generateData(false);  % false = no calcular SVD

% Exportar datos estándar
outputDir = 'Tutorials/NilHeinzFiles';
fprintf('\nExportando datos estándar a CSV...\n');
dataGen.exportToCSV(outputDir);
fprintf('  ✓ TTrainingData.csv\n');
fprintf('  ✓ KcoarseTrainingData.csv\n');
fprintf('  ✓ McoarseTrainingData.csv\n');

% Verificar que SVD NO está calculado
fprintf('\nVerificación: SVD calculado? %s\n', ...
    mat2str(~isempty(dataGen.getSVDResults())));

%% ================================================================
%% FLUJO 2: Calcular SVD de forma independiente DESPUÉS
%% ================================================================
fprintf('\n--- FLUJO 2: Cálculo independiente de SVD ---\n');

% Calcular SVD a partir de los datos ya generados
fprintf('\nCalculando SVD a partir de datos existentes...\n');
dataGen.computeSVDFromExistingData();

% Verificar que SVD ahora SÍ está calculado
[U, S, V, k] = dataGen.getSVDResults();
fprintf('\nVerificación: SVD calculado? %s\n', ...
    mat2str(~isempty(U)));
fprintf('  - Número de modos retenidos: %d\n', k);
fprintf('  - Dimensiones U: [%d × %d]\n', size(U,1), size(U,2));
fprintf('  - Dimensiones S: [%d × %d]\n', size(S,1), size(S,2));
fprintf('  - Dimensiones V: [%d × %d]\n', size(V,1), size(V,2));

% Exportar datos SVD
fprintf('\nExportando datos SVD...\n');
dataGen.exportSVDToCSV(outputDir);
dataGen.exportSVDToMAT(outputDir, 'SVD_Results_Independent.mat');
fprintf('  ✓ VTrainingData.csv\n');
fprintf('  ✓ SVD_Results_Independent.mat\n');

%% ================================================================
%% FLUJO 3: Verificar que se puede recalcular SVD con diferentes parámetros
%% ================================================================
fprintf('\n--- FLUJO 3: Recalcular SVD con diferentes parámetros ---\n');

% Crear nuevo generador con parámetros SVD más estrictos
dataGen2 = TrainingDataGenerator(radii, ...
    'svdTolerance', 1e-8, ...      % Más estricto
    'svdEnergyRatio', 0.999);      % Más energía requerida

% Copiar datos generados (simulando que ya los tenemos)
dataGen2.allResults = dataGen.allResults;
dataGen2.TData = dataGen.TData;
dataGen2.KData = dataGen.KData;
dataGen2.MData = dataGen.MData;

% Recalcular SVD con nuevos parámetros
fprintf('\nRecalculando SVD con parámetros más estrictos...\n');
dataGen2.computeSVDFromExistingData();

k2 = dataGen2.getNumberOfModes();
fprintf('  - Modos retenidos (tolerancia 1e-6): %d\n', k);
fprintf('  - Modos retenidos (tolerancia 1e-8): %d\n', k2);

%% ================================================================
%% FLUJO 4: Comparación directa - Todo en un paso
%% ================================================================
fprintf('\n--- FLUJO 4: Comparación - Todo en un paso ---\n');

% Crear nuevo generador y generar TODO (datos + SVD) en un paso
dataGen3 = TrainingDataGenerator(radii);
fprintf('\nGenerando datos Y SVD en un solo paso...\n');
dataGen3.generateData(true);  % true = calcular SVD también

k3 = dataGen3.getNumberOfModes();
fprintf('  - Modos retenidos: %d\n', k3);

% Verificar que ambos métodos dan el mismo resultado
if k == k3
    fprintf('\n✓ Ambos métodos (independiente vs. integrado) dan el mismo resultado.\n');
else
    fprintf('\n⚠ Diferencia en número de modos (puede ser por parámetros diferentes).\n');
end

%% ================================================================
%% Resumen
%% ================================================================
fprintf('\n=== RESUMEN ===\n');
fprintf('✓ Flujo 1: Generación sin SVD - COMPLETADO\n');
fprintf('✓ Flujo 2: Cálculo independiente de SVD - COMPLETADO\n');
fprintf('✓ Flujo 3: Recalcular SVD con diferentes parámetros - COMPLETADO\n');
fprintf('✓ Flujo 4: Generación integrada (datos + SVD) - COMPLETADO\n');
fprintf('\nLa clase permite ambos enfoques:\n');
fprintf('  1. generateData(false) → computeSVDFromExistingData()\n');
fprintf('  2. generateData(true)  (todo en un paso)\n');
fprintf('\nAmbos son equivalentes y permiten flexibilidad en el flujo de trabajo.\n');

