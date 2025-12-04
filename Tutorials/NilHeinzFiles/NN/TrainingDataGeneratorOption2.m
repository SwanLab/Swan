% EXAMPLE_OPTIMIZEDDATAGENERATION Ejemplo de uso de OptimizedDataGeneration
 close all;

%% Configuración
radii = 0.5*(0.1:0.1:0.9);  % Radius to include as training set

%% Crear generador de datos
dataGen = TrainingDataGenerator(radii);

%% Generar datos
dataGen.generateData(false); 

%% Exportar
outputDir = [userpath '\Tutorials\NilHeinzFiles\NN\TrainingDataFiles'];
dataGen.exportToCSV(outputDir);
   