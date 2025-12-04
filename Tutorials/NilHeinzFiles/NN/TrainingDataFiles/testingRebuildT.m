
% Crear generador de datos

load('referenceMesh.mat');

[K, M, T] = RebuildKMTData.compute('A', 0.25, mR);