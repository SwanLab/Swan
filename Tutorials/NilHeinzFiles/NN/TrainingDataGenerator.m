classdef TrainingDataGenerator < handle
    % TRAININGDATAGENERATOR Generación optimizada de datos para entrenamiento de NNs
    %
    % Esta clase reutiliza Training.m y OfflineDataProcessor.m para generar
    % datos de entrenamiento de forma eficiente.
    
    properties (Access = private)
        % Parámetros de generación
        radii              % Vector de radios a procesar
               
        % Datos generados
        TData              % Datos de T: [r, x, y, Tx1, Ty1, ..., Tx8, Ty8]
        KData              % Datos de K: [r, K11, K12, ..., K88]
        MData              % Datos de M: [r, M11, M12, ..., M88]
        
        % Resultados intermedios
        allResults         % Cell array con resultados por radio
        
        % Datos SVD
        T_SVD              % Matriz T_SVD: [nDofs×8 × nRadii]
        U                  % Modos espaciales: [nDofs×8 × k]
        S                  % Valores singulares: [k × k] (diagonal)
        V                  % Modos paramétricos: [nRadii × k]
        k                  % Número de modos retenidos
        VTrainingData      % Datos de entrenamiento para V: [r, V1, V2, ..., Vk]
        svdTolerance       % Tolerancia para truncamiento SVD (default: 1e-6)
        svdEnergyRatio     % Ratio de energía para truncamiento (default: 0.99)
    end
    
    methods (Access = public)
        
        function obj = TrainingDataGenerator(radii, varargin)
            % Constructor
            % radii: Vector de radios (ej: 0:0.1:0.9)
            % varargin: 
            %   'svdTolerance': Tolerancia para truncamiento SVD (default: 1e-6)
            %   'svdEnergyRatio': Ratio de energía para truncamiento (default: 0.99)
            
            obj.radii = radii;
            
            % Parsear parámetros opcionales
            p = inputParser;
            addParameter(p, 'svdTolerance', 1e-6, @isnumeric);
            addParameter(p, 'svdEnergyRatio', 0.99, @isnumeric);
            parse(p, varargin{:});
            
            obj.svdTolerance = p.Results.svdTolerance;
            obj.svdEnergyRatio = p.Results.svdEnergyRatio;
        end
        
        function generateData(obj, computeSVD)
            % Genera todos los datos de entrenamiento
            % computeSVD: (opcional) true para calcular SVD, false para omitir (default: false)
            %
            % Nota: Si los datos ya están generados, solo calculará SVD si computeSVD=true
            %       y SVD aún no ha sido calculado. Para recalcular SVD, use computeSVDFromExistingData()
            
            if nargin < 2
                computeSVD = false;
            end
            
            % Verificar si los datos ya están generados
            if ~isempty(obj.allResults) && length(obj.allResults) == length(obj.radii)
                fprintf('Datos ya generados. Saltando generación...\n');
            else
                fprintf('Generando datos para %d radios...\n', length(obj.radii));
                
                % Pre-allocar resultados
                nRadii = length(obj.radii);
                obj.allResults = cell(nRadii, 1);
                
                for j = 1:nRadii
                    obj.allResults{j} = obj.solveForRadius(obj.radii(j), j);
                    if mod(j, 10) == 0
                        fprintf('  Procesado %d/%d radios\n', j, nRadii);
                    end
                end
                
                % Post-procesamiento
                obj.processTData();
                obj.processKMData();
            end
            
            % Calcular SVD si se solicita y aún no está calculado
            if computeSVD && isempty(obj.U)
                fprintf('Calculando SVD...\n');
                obj.computeSVD();
                obj.processSVDData();
                fprintf('SVD completado. Modos retenidos: %d\n', obj.k);
            elseif computeSVD && ~isempty(obj.U)
                fprintf('SVD ya calculado. Use computeSVDFromExistingData() para recalcular.\n');
            end
            
            fprintf('Generación completada.\n');
        end
        
        function exportToCSV(obj, outputDir)
            % Exporta datos a archivos CSV
            % outputDir: Directorio de salida
            
            % Export datos estándar
            TFile = fullfile(outputDir, 'TTrainingData.csv');
            writematrix(obj.TData, TFile);

            McFile = fullfile(outputDir, 'McoarseTrainingData.csv');
            writematrix(obj.MData, McFile);
            
            KcFile = fullfile(outputDir, 'KcoarseTrainingData.csv');
            writematrix(obj.KData, KcFile);
        end
        
        function exportSVDToCSV(obj, outputDir)
            % Exporta datos de entrenamiento SVD a CSV
            % outputDir: Directorio de salida
            % Genera: VTrainingData.csv con formato [r, V1, V2, ..., Vk]
            
            if isempty(obj.VTrainingData)
                error('SVD no ha sido calculado. Llame a generateData(true) o computeSVDFromExistingData() primero.');
            end
            
            VFile = fullfile(outputDir, 'VTrainingData.csv');
            writematrix(obj.VTrainingData, VFile);
            fprintf('Datos SVD exportados a: %s\n', VFile);
        end
        
        function exportSVDToMAT(obj, outputDir, fileName)
            % Exporta resultados SVD a archivo .mat para reconstrucción
            % outputDir: Directorio de salida
            % fileName: Nombre del archivo (default: 'SVD_Results.mat')
            
            if isempty(obj.U) || isempty(obj.S) || isempty(obj.V)
                error('SVD no ha sido calculado. Llame a generateData(true) o computeSVDFromExistingData() primero.');
            end
            
            if nargin < 3
                fileName = 'SVD_Results.mat';
            end
            
            % Obtener mesh de referencia (del primer resultado)
            if ~isempty(obj.allResults)
                meshRef = obj.allResults{1}.mesh;
            else
                meshRef = [];
            end
            
            % Guardar resultados SVD
            U = obj.U;
            S = obj.S;
            V = obj.V;
            k = obj.k;
            radii = obj.radii;
            
            filePath = fullfile(outputDir, fileName);
            save(filePath, 'U', 'S', 'V', 'k', 'radii', 'meshRef', '-v7.3');
            fprintf('Resultados SVD exportados a: %s\n', filePath);
            fprintf('  - U: [%d × %d] modos espaciales\n', size(U,1), size(U,2));
            fprintf('  - S: [%d × %d] valores singulares\n', size(S,1), size(S,2));
            fprintf('  - V: [%d × %d] modos paramétricos\n', size(V,1), size(V,2));
            fprintf('  - k: %d modos retenidos\n', k);
        end
        
        function k = getNumberOfModes(obj)
            % Retorna el número de modos retenidos
            k = obj.k;
        end
        
        function [U, S, V, k] = getSVDResults(obj)
            % Retorna los resultados SVD
            U = obj.U;
            S = obj.S;
            V = obj.V;
            k = obj.k;
        end
        
        function computeSVDFromExistingData(obj)
            % Calcula SVD a partir de datos ya generados
            % Útil cuando se generaron datos sin SVD y luego se quiere calcular SVD
            % independientemente
            
            if isempty(obj.allResults)
                error('No hay datos generados. Llame a generateData() primero.');
            end
            
            fprintf('Calculando SVD a partir de datos existentes...\n');
            obj.computeSVD();
            obj.processSVDData();
            fprintf('SVD completado. Modos retenidos: %d\n', obj.k);
        end
         
    end
    
    methods (Access = private)
        
        function mesh = createMesh(obj,r)
            
            fullmesh = UnitTriangleMesh(12,12);
            ls = obj.computeCircleLevelSet(fullmesh,r);
            sUm.backgroundMesh = fullmesh;
            sUm.boundaryMesh   = fullmesh.createBoundaryMesh;
            uMesh              = UnfittedMesh(sUm);
            uMesh.compute(ls);
            holeMesh = uMesh.createInnerMesh();
            s.coord     = holeMesh.coord;
            s.connec    = holeMesh.connec;
            s.interType = 'LINEAR';
            s           = obj.updateCoordsMesh(s);
            mesh = Mesh.create(s);

        end
        
        function ls = computeCircleLevelSet(obj, mesh, radius)
            % Original method for single radius (backward compatibility)
            gPar.type          = 'Circle';
            gPar.radius        = radius;
            gPar.xCoorCenter   = 0.5;
            gPar.yCoorCenter   = 0.5;
            g                  = GeometricalFunction(gPar);
            phiFun             = g.computeLevelSetFunction(mesh);
            lsCircle           = phiFun.fValues;
            ls = -lsCircle;
        end

        function s = updateCoordsMesh(obj, s)
            % Nudge nodes at the four rectangle corners in x to avoid
            % exact coincidences.
            tol  = 1e-8;
            epsx = 1e-9;
        
            x = s.coord(:,1); y = s.coord(:,2);
            xmax = max(x); xmin = min(x);
            ymax = max(y); ymin = min(y);
        
            % Top-right (xmax,ymax)
            mask = abs(x - xmax) < tol & abs(y - ymax) < tol;
            s.coord(mask, :) = s.coord(mask, :) - [epsx, 0];
        
            % Bottom-right (xmax,ymin)
            mask = abs(x - xmax) < tol & abs(y - ymin) < tol;
            s.coord(mask, :) = s.coord(mask, :) - [epsx, 0];
        
            % Top-left (xmin,ymax)
            mask = abs(x - xmin) < tol & abs(y - ymax) < tol;
            s.coord(mask, :) = s.coord(mask, :) + [epsx, 0];
        
            % Bottom-left (xmin,ymin)
            mask = abs(x - xmin) < tol & abs(y - ymin) < tol;
            s.coord(mask, :) = s.coord(mask, :) + [epsx, 0];
        end
        
        function result = solveForRadius(obj, r, idx)
            % Resuelve el problema elástico para un radio dado usando Training
            % r: Radio de la inclusión
            % idx: Índice del radio (para logging)
            
            % Crear malla base
            meshRef = obj.createMesh(r);
                 
            trainingData = Training(meshRef);  % Con radio, material variable
            
            % Extraer datos: uSbd tiene 8 columnas (una por modo), LHSsbd es K_fine
            u = trainingData.uSbd;  % [nDofs × 8] - Esta es la matriz T
            
            % Calcular Kcoarse usando OfflineDataProcessor
            processor = OfflineDataProcessor(trainingData);
            EIFEoper = processor.computeROMbasis(r);
            Kcoarse = EIFEoper.Kcoarse;
            Mcoarse = EIFEoper.Mcoarse;
            
            % Almacenar resultado
            result.r = r;
            result.u = u;
            result.Kcoarse = Kcoarse;
            result.Mcoarse = Mcoarse;
            result.mesh = trainingData.mesh;
        end
        
        function processTData(obj)
            % Procesa y formatea datos de T para CSV
            % Formato: [r, x, y, Tx1, Ty1, Tx2, Ty2, ..., Tx8, Ty8]
            
            TData = [];
            
            for j = 1:length(obj.allResults)
                result = obj.allResults{j};
                r = result.r;
                u = result.u;  % [nDofs × 8]
                mesh = result.mesh;
                
                % Reshape u: cada columna tiene [u_x1, u_y1, u_x2, u_y2, ...]
                % Necesitamos [Tx1, Ty1, Tx2, Ty2, ...] por nodo
                nnodes = mesh.nnodes;
                ndim = mesh.ndim;
                
                % Para cada modo, extraer componentes x e y
                t_reshaped = zeros(nnodes, 16);  % 8 modos × 2 componentes
                for mode = 1:8
                    u_mode = u(:, mode);  % [nDofs × 1]
                    u_reshaped = reshape(u_mode, ndim, nnodes)';  % [nnodes × 2]
                    t_reshaped(:, 2*mode-1:2*mode) = u_reshaped;  % [Tx_mode, Ty_mode]
                end
                
                % Combinar: [r, x, y, Tx1, Ty1, ..., Tx8, Ty8]
                t_aux = [r * ones(nnodes, 1), mesh.coord, t_reshaped];
                TData = [TData; t_aux];
            end
            
            obj.TData = TData;
        end
        
        function processKMData(obj)
            % Procesa y formatea datos de K para CSV
            % Formato: [r, K11, K12, K13, ..., K88] (36 componentes únicas)
            
            nRadii = length(obj.radii);
            Kdata = zeros(nRadii, 36);  % 36 componentes únicas de K 8×8
            Mdata = zeros(nRadii, 36);
            for j = 1:nRadii
                result = obj.allResults{j};
                K = result.Kcoarse;  % [8 × 8]
                M = result.Mcoarse;
                % Extraer triangular superior
                triangSupK = triu(K);
                triangSupM = triu(M);
                rowK = [];
                rowM = [];
                for i = 1:8
                    for k = i:8
                        rowK(end+1) = triangSupK(i, k);
                        rowM(end+1) = triangSupM(i, k);
                    end
                end
                Kdata(j, :) = rowK;
                Mdata(j, :) = rowM;
            end
            
            % Agregar columna de radios
            obj.KData = [obj.radii(:), Kdata];
            obj.MData = [obj.radii(:), Mdata];
        end
        
        function computeSVD(obj)
            % Construye T_SVD y aplica descomposición SVD
            % Reutiliza los datos T ya generados en allResults
            
            nRadii = length(obj.allResults);
            
            % Verificar que todos los resultados tienen la misma dimensión
            firstResult = obj.allResults{1};
            u_first = firstResult.u;  % [nDofs × 8]
            nDofs = size(u_first, 1);
            
            % Pre-allocar T_SVD: [nDofs×8 filas × nRadii columnas]
            obj.T_SVD = zeros(nDofs * 8, nRadii);
            
            % Construir T_SVD: cada columna es T vectorizado para un radio
            for j = 1:nRadii
                result = obj.allResults{j};
                u = result.u;  % [nDofs × 8]
                % Vectorizar: convertir [nDofs × 8] → [nDofs×8 × 1]
                obj.T_SVD(:, j) = u(:);
            end
            
            fprintf('  T_SVD construido: [%d × %d]\n', size(obj.T_SVD, 1), size(obj.T_SVD, 2));
            
            % Aplicar SVD
            fprintf('  Aplicando SVD...\n');
            [U_full, S_full, V_full] = svd(obj.T_SVD, 'econ');
            
            % Determinar número de modos a retener
            obj.k = obj.determineNumberOfModes(S_full);
            
            % Truncar a k modos
            obj.U = U_full(:, 1:obj.k);
            obj.S = S_full(1:obj.k, 1:obj.k);
            obj.V = V_full(:, 1:obj.k);
            
            fprintf('  SVD completado. Modos retenidos: %d/%d\n', obj.k, min(size(obj.T_SVD)));
        end
        
        function k = determineNumberOfModes(obj, S)
            % Determina el número de modos a retener basado en tolerancia y energía
            % S: Matriz de valores singulares (diagonal)
            
            sigma = diag(S);
            sigma1 = sigma(1);
            
            % Criterio 1: Tolerancia relativa
            k_tol = find(sigma / sigma1 > obj.svdTolerance, 1, 'last');
            if isempty(k_tol)
                k_tol = 1;
            end
            
            % Criterio 2: Ratio de energía
            energy = cumsum(sigma.^2) / sum(sigma.^2);
            k_energy = find(energy >= obj.svdEnergyRatio, 1);
            if isempty(k_energy)
                k_energy = length(sigma);
            end
            
            % Usar el mínimo de ambos criterios (más conservador)
            k = min(k_tol, k_energy);
            
            % Asegurar al menos 1 modo y no más que el máximo disponible
            k = max(1, min(k, length(sigma)));
            
            % Información de diagnóstico
            fprintf('    Criterio tolerancia: %d modos (σ_%d/σ_1 = %.2e)\n', ...
                k_tol, k_tol, sigma(k_tol)/sigma1);
            fprintf('    Criterio energía: %d modos (%.1f%% energía)\n', ...
                k_energy, energy(k_energy)*100);
            fprintf('    Modos seleccionados: %d\n', k);
        end
        
        function processSVDData(obj)
            % Procesa datos SVD para entrenamiento de redes neuronales
            % Genera VTrainingData: [r, V1, V2, ..., Vk]
            % Cada fila corresponde a un radio, cada columna a un modo paramétrico
            
            if isempty(obj.V)
                error('SVD no ha sido calculado. Llame a computeSVD() primero.');
            end
            
            % V ya tiene el formato correcto: [nRadii × k]
            % Agregar columna de radios
            obj.VTrainingData = [obj.radii(:), obj.V];
            
            fprintf('  Datos de entrenamiento V generados: [%d × %d]\n', ...
                size(obj.VTrainingData, 1), size(obj.VTrainingData, 2));
        end
        
        function T_reconstructed = reconstructT(obj, r_index)
            % Reconstruye T para un radio dado usando SVD truncado
            % r_index: Índice del radio en obj.radii (o puede ser el radio mismo)
            
            if isempty(obj.U) || isempty(obj.S) || isempty(obj.V)
                error('SVD no ha sido calculado. Llame a computeSVD() primero.');
            end
            
            % Si r_index es un radio, encontrar su índice
            if isscalar(r_index) && r_index > 0 && r_index <= 1
                [~, idx] = min(abs(obj.radii - r_index));
                r_index = idx;
            end
            
            if r_index < 1 || r_index > size(obj.V, 1)
                error('Índice de radio fuera de rango.');
            end
            
            % Reconstrucción: T = U × S × V'
            V_row = obj.V(r_index, :)';  % [k × 1]
            T_vec = obj.U * obj.S * V_row;  % [nDofs×8 × 1]
            
            % Reshape a [nDofs × 8]
            firstResult = obj.allResults{1};
            u_first = firstResult.u;
            nDofs = size(u_first, 1);
            T_reconstructed = reshape(T_vec, [nDofs, 8]);
        end
        
        function error = computeReconstructionError(obj, r_index)
            % Calcula el error de reconstrucción SVD para un radio
            % error: Error relativo ||T_original - T_reconstructed|| / ||T_original||
            
            if r_index < 1 || r_index > length(obj.allResults)
                error('Índice de radio fuera de rango.');
            end
            
            % T original
            T_original = obj.allResults{r_index}.u;  % [nDofs × 8]
            
            % T reconstruido
            T_reconstructed = obj.reconstructT(r_index);
            
            % Error relativo
            error = norm(T_original - T_reconstructed, 'fro') / norm(T_original, 'fro');
        end
        
    end
    
end
