% --- Script para extraer datos de subplots específicos de varios .fig ---
clear all, close all, clc;

% Lista de archivos .fig que quieres procesar
files = {'0_3D_LevelSet_try3.fig','Monitoring_45_3D_LevelSet.fig', 'Monitoring_90_3D_LevelSet.fig','Monitoring_0_90_3D_LevelSet',...
         'Monitoring_45_45_3D_LevelSet.fig', 'Monitoring_0_45_3D_LevelSet.fig'};

% Inicializamos estructura para guardar los datos
allData = struct();

for f = 1:length(files)
    % Abrimos el .fig (sin mostrar)
    fig = openfig(files{f}, 'invisible');
    
    % Obtenemos todos los ejes (subplots) en orden
    allAxes = findall(fig, 'type', 'axes');
    
    % IMPORTANTE: en MATLAB los axes suelen estar en orden inverso
    % Así que los reordenamos por posición
    allAxes = flipud(allAxes);
    
    % Grid: 2 filas × 5 columnas
    nRows = 2; nCols = 5;
    
    % Queremos fila 1, columnas 2 y 3
    targetSubs = [ sub2ind([nCols, nRows], 2, 1), ...
                   sub2ind([nCols, nRows], 3, 1) ];
    % Nota: sub2ind usa formato (nCols,nRows), ojo al orden
    
    % Inicializamos
    figData = {};
    
    for t = 1:length(targetSubs)
        ax = allAxes(targetSubs(t));
        lines = findall(ax, 'Type', 'Line');
        
        % En caso de que haya varias líneas, cogemos la primera
        if ~isempty(lines)
            x = get(lines(1), 'XData');
            y = get(lines(1), 'YData');
            figData{end+1} = struct('x', x, 'y', y);
        end
    end
    
    % Guardamos en estructura con nombre de archivo
    allData.(sprintf('fig%d', f)) = figData;
    
    close(fig); % cerramos el archivo
end

%% --- Ejemplo de replotear los vectores de todos los .fig ---
% Subplot 1 (columna 2, fila 1)
figure(1); hold on; grid minor;
for f = 1:length(files)
    data = allData.(sprintf('fig%d', f)){1}; % primer subplot
    plot(data.x, data.y, 'LineWidth',1.5);
end
xlabel('Iteration')
%ylim([1, 3])
legend('0º', '45º','90º','0º/90º','-45º/45º','0º/45º');
title('Compliance comparison - Cantilever, LevelSet')

% Subplot 2 (columna 3, fila 1)
figure(2); hold on; grid minor;
for f = 1:length(files)
    data = allData.(sprintf('fig%d', f)){2}; % segundo subplot
    plot(data.x, data.y, 'LineWidth',1.5);
end
xlabel('Iteration')
legend('0º', '45º','90º','0º/90º','45º/45º','0º/45º');
title('Volume Constraint comparison - Cantilever, LevelSet')