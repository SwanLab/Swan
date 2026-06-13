%% Script corregido: Malla completa vs. Celdas de referencia
clear; clc; close all;

% Dimensiones del subdominio base
numSubX = 4; 
numSubY = 2;
globalLength = 2.0;
globalHeight = 0.5;
subLength = globalLength / numSubX;
subHeight = globalHeight / numSubY;

%% 1. Cargar y procesar los datos EXACTAMENTE como tus funciones

% --- PANEL 1: Malla Estructurada Completa (5x5 elementos por subdominio) ---
nPerSide = 6; 
V_struct = []; F_struct = []; node_offset = 0;
for i = 1:numSubX
    for j = 1:numSubY
        x_start = (i-1) * subLength;
        y_start = (j-1) * subHeight;
        x1 = linspace(x_start, x_start + subLength, nPerSide);
        x2 = linspace(y_start, y_start + subHeight, nPerSide);
        [xv, yv] = meshgrid(x1, x2);
        [F, V] = mesh2tri(xv, yv, zeros(size(xv)), 'x');
        F_struct = [F_struct; F + node_offset];
        V_struct = [V_struct; V(:, 1:2)];
        node_offset = size(V_struct, 1);
    end
end

% --- PANEL 2: Un SOLO elemento Lattice de referencia ---
data_lat = load('mallaLattice2D.mat');
campos_lat = fieldnames(data_lat);
meshData_lat = data_lat.(campos_lat{1});
V_lattice = meshData_lat.coord;
F_lattice = meshData_lat.connec;
if size(F_lattice, 2) == 4
    F_lattice = [F_lattice(:, [1 2 3]); F_lattice(:, [1 3 4])];
end
% Escalado idéntico a tu función createLatticeMesh
minX = min(V_lattice(:,1)); maxX = max(V_lattice(:,1));
minY = min(V_lattice(:,2)); maxY = max(V_lattice(:,2));
escalaGlobal = min(subLength / (maxX - minX), subHeight / (maxY - minY));
V_lattice(:,1) = (V_lattice(:,1) - minX) * escalaGlobal;
V_lattice(:,2) = (V_lattice(:,2) - minY) * escalaGlobal;

% --- PANEL 3: Un SOLO elemento Auxético de referencia ---
data_aux = load('DEF_Q4auxL_1.mat');
V_aux = data_aux.EIFEoper.MESH.COOR;
cnQ4 = double(data_aux.EIFEoper.MESH.CN);
F_aux = [cnQ4(:, [1 2 3]); cnQ4(:, [1 3 4])];
% Escalado idéntico a tu función createAuxeticMesh
minX = min(V_aux(:,1)); maxX = max(V_aux(:,1));
minY = min(V_aux(:,2));
scale_aux = subLength / (maxX - minX);
V_aux(:,1) = (V_aux(:,1) - minX) * scale_aux;
V_aux(:,2) = (V_aux(:,2) - minY) * scale_aux;


%% 2. Graficado con Tiledlayout (Estructurada arriba, referencias abajo)

figure('Color', 'w', 'Units', 'inches', 'Position', [1, 1, 7, 5]);
t = tiledlayout(2, 2, 'TileSpacing', 'Loose', 'Padding', 'Compact');

% --- Fila superior: Malla Estructurada completa (Ocupa las 2 columnas) ---
nexttile(t, [1, 2]);
hold on; box on;
patch('Faces', F_struct, 'Vertices', V_struct, ...
      'FaceColor', [0.97 0.97 0.97], 'EdgeColor', [0.6 0.6 0.6], 'LineWidth', 0.4);
% Líneas de subdominios internas
for i = 1:numSubX-1
    plot([i*subLength, i*subLength], [0, globalHeight], 'k-', 'LineWidth', 1.3);
end
for j = 1:numSubY-1
    plot([0, globalLength], [j*subHeight, j*subHeight], 'k-', 'LineWidth', 1.3);
end
plot([0, globalLength, globalLength, 0, 0], [0, 0, globalHeight, globalHeight, 0], 'k-', 'LineWidth', 1.6);
title('(a) Global Structured Mesh (5x5 elements per subdomain)', 'FontSize', 10.5, 'FontWeight', 'bold');
axis equal; 
xlim([-0.05, globalLength + 0.05]); ylim([-0.05, globalHeight + 0.05]);
set(gca, 'XTick', [0, globalLength], 'YTick', [0, globalHeight], 'FontSize', 9);

% --- Fila inferior, Izquierda: Celda Unidad Lattice ---
nexttile(t);
hold on; box on;
patch('Faces', F_lattice, 'Vertices', V_lattice, ...
      'FaceColor', [0.95 0.95 0.95], 'EdgeColor', [0.5 0.5 0.5], 'LineWidth', 0.5);
plot([0, subLength, subLength, 0, 0], [0, 0, subHeight, subHeight, 0], 'k--', 'LineWidth', 1.2);
title('(b) Lattice Reference Cell', 'FontSize', 10.5, 'FontWeight', 'bold');
axis equal;
xlim([-0.02, subLength + 0.02]); ylim([-0.02, subHeight + 0.02]);
xlabel(['l = ', num2str(subLength)], 'FontSize', 9); ylabel(['h = ', num2str(subHeight)], 'FontSize', 9);
set(gca, 'XTick', [0, subLength], 'YTick', [0, subHeight], 'FontSize', 8);

% --- Fila inferior, Derecha: Celda Unidad Auxética ---
nexttile(t);
hold on; box on;
patch('Faces', F_aux, 'Vertices', V_aux, ...
      'FaceColor', [0.95 0.95 0.95], 'EdgeColor', [0.5 0.5 0.5], 'LineWidth', 0.5);
plot([0, subLength, subLength, 0, 0], [0, 0, subHeight, subHeight, 0], 'k--', 'LineWidth', 1.2);
title('(c) Auxetic Reference Cell', 'FontSize', 10.5, 'FontWeight', 'bold');
axis equal;
xlim([-0.02, subLength + 0.02]); ylim([-0.02, subHeight + 0.02]);
xlabel(['l = ', num2str(subLength)], 'FontSize', 9);
set(gca, 'XTick', [0, subLength], 'YTick', [0, subHeight], 'FontSize', 8);

%% 3. Exportar
exportgraphics(gcf, 'reference_meshes_2d.pdf', 'ContentType', 'vector');