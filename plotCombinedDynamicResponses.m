% plotCombinedDynamicResponses.m
% Script para combinar varias gráficas de vibraciones en plots únicos

clc; close all;

% =========================================================================
% CONFIGURACIÓN DE CASOS
% =========================================================================
% Definimos los casos a procesar: Nombre, Título, y lista de archivos
cases = {
    'STEP', 'Combined Dynamic Response (STEP)', {
        'STEPMAT1.fig', 'Material 1 (IP)';
        'STEPMAT2.fig', 'Material 2 (SL)';
        'STEPMAT4.fig', 'Material 4 (MDL)'
    };
    'SINUSOIDAL', 'Combined Dynamic Response (SINUSOIDAL)', {
        'SINMAT1.fig', 'Material 1 (IP)';
        'SINMAT2.fig', 'Material 2 (SL)';
        'SINMAT4.fig', 'Material 4 (MDL)'
    }
};

% =========================================================================
% GENERACIÓN DE LOS PLOTS COMBINADOS
% =========================================================================
% Resetear ajustes de texto (por si acaso)
set(groot, 'defaultTextInterpreter', 'remove');
set(groot, 'defaultAxesTickLabelInterpreter', 'remove');
set(groot, 'defaultLegendInterpreter', 'remove');

% Extraer los colores por defecto de MATLAB (parula palette)
colors = colororder; 

for c = 1:size(cases, 1)
    caseName = cases{c, 1};
    titleStr = cases{c, 2};
    figFiles = cases{c, 3};
    
    % Crear nueva figura para el caso actual
    hCombined = figure('Name', sprintf('Combined Dynamic Response - %s', caseName));
    hAx = axes(hCombined);
    hold(hAx, 'on');
    grid(hAx, 'on'); box(hAx, 'on');
    
    for i = 1:size(figFiles, 1)
        filename = figFiles{i, 1};
        labelStr = figFiles{i, 2};
        
        if ~isfile(filename)
            warning('No se ha encontrado el archivo %s. Saltando...', filename);
            continue;
        end
        
        % Abrir figura de forma invisible
        figTemp = openfig(filename, 'invisible');
        
        % Buscar las líneas dentro de la figura
        lineObjs = findobj(figTemp, 'Type', 'Line');
        
        % Identificar la línea con los datos reales (la del 0 solo tiene 2 puntos)
        dataLine = [];
        for j = 1:length(lineObjs)
            if length(lineObjs(j).XData) > 2
                dataLine = lineObjs(j);
                break;
            end
        end
        
        if ~isempty(dataLine)
            % Asignar el GRANATE por defecto de MATLAB (el 7º color) si es el Material 4
            if contains(labelStr, 'Material 4')
                thisColor = colors(7, :);
            else
                thisColor = colors(i, :);
            end
            
            % Volver a plotear usando los datos extraídos (convertidos a mm)
            plot(hAx, dataLine.XData, dataLine.YData * 1000, '-', 'LineWidth', 2.0, ...
                'Color', thisColor, 'DisplayName', labelStr);
        else
            warning('No se pudo extraer la curva de %s', filename);
        end
        
        % Cerrar la figura temporal
        close(figTemp);
    end
    
    % Añadir la línea horizontal negra de y=0
    xlims = hAx.XLim;
    plot(hAx, xlims, [0, 0], 'k--', 'LineWidth', 1.0, 'HandleVisibility', 'off');
    
    % Aplicar el formato estético LaTeX
    xlabel(hAx, 'Time $t$ (s)', 'FontSize', 13, 'Interpreter', 'latex');
    ylabel(hAx, 'Displacement $w$ (mm)', 'FontSize', 13, 'Interpreter', 'latex');
    title(hAx, titleStr, 'FontSize', 14, 'Interpreter', 'latex');
    set(hAx, 'FontSize', 13, 'GridAlpha', 0.3, 'TickLabelInterpreter', 'latex');
    
    % Añadir leyenda interactiva
    legend(hAx, 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 12);
    hold(hAx, 'off');
end

fprintf('¡Plots combinados (STEP y SINUSOIDAL) generados exitosamente!\n');
