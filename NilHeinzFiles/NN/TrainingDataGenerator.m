classdef TrainingDataGenerator < handle
    % TRAININGDATAGENERATOR Generación optimizada de datos para entrenamiento de NNs
    
    properties (Access = private)
        % Parámetros de generación
        radii              
        referenceMesh      % Mesh de referencia fijo para todos los radios
        referenceRadius    % Radio usado para crear el mesh de referencia
               
        % Datos generados
        TData              % Datos de T: [r, x, y, Tx1, Ty1, ..., Tx8, Ty8]
        KData              % Datos de K: [r, K11, K12, ..., K88]
        MData              % Datos de M: [r, M11, M12, ..., M88]
        
        % Resultados intermedios
        allResults         % Cell array con resultados por radio
        
        % Datos SVD
        T_SVD              % Matriz T_SVD: [nDofs×8 × nRadii]
        U_full             % FUll Spatial Modes
        U                  % Spatial Modes Keeping K modes [nDofs×8 × k]
        S                  % Valores singulares: [k × k] (diagonal)
        V                  % Modos paramétricos: [nRadii × k]
        k                  % Número de modos retenidos
        VTrainingData      % Datos de entrenamiento para V: [r, V1, V2, ..., Vk]
    end
    
    methods (Access = public)
        
        function obj = TrainingDataGenerator(radii, varargin)

            obj.radii = radii;
            
            % Parsear parámetros opcionales
            p = inputParser;
            addParameter(p, 'referenceRadius', max(radii), @isnumeric);
            parse(p, varargin{:});
            
            obj.referenceRadius = p.Results.referenceRadius;
            
            % Crear mesh de referencia una vez para garantizar dimensiones consistentes
            obj.referenceMesh = obj.createMesh(obj.referenceRadius);
        end
        
        function generateData(obj, computeSVD)
            % Genera todos los datos de entrenamiento
           
            fprintf('Generating training data for %d radius...\n', length(obj.radii));
                
                % Pre-allocar resultados
                nRadii = length(obj.radii);
                obj.allResults = cell(nRadii, 1);
                
                for j = 1:nRadii
                    obj.allResults{j} = obj.solveForRadius(obj.radii(j));
                    if mod(j, 10) == 0
                        fprintf('  Procesado %d/%d radios\n', j, nRadii);
                    end
                end
                
                % Post-procesamiento
                obj.processTData();
                obj.processKMData();
            
            % SVD if requested
            if computeSVD
                fprintf('Calculando SVD...\n');
                obj.computeSVD();
                obj.processSVDData();
                fprintf('SVD completado. Modos retenidos: %d\n', obj.k);
            end
            
            fprintf('Generación completada.\n');
        end
        
        function exportToCSV(obj, outputDir)
            % Exporta datos a archivos CSV
            
            TFile = fullfile(outputDir, 'TTrainingData.csv');
            writematrix(obj.TData, TFile);

            McFile = fullfile(outputDir, 'McoarseTrainingData.csv');
            writematrix(obj.MData, McFile);
            
            KcFile = fullfile(outputDir, 'KcoarseTrainingData.csv');
            writematrix(obj.KData, KcFile);
        end
        
        function exportSVDToCSV(obj, outputDir)
            % Genera: VTrainingData.csv con formato [r, V1, V2, ..., Vk]
            
            VFile = fullfile(outputDir, 'VTrainingData.csv');
            writematrix(obj.VTrainingData, VFile);
            fprintf('Datos SVD exportados a: %s\n', VFile);
        end
        
        function exportSVDToMAT(obj,fileName, outputDir)
            
            
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
            k = obj.k;
        end        
         
    end
    
    methods (Access = private)
        
        function mesh = createMesh(obj,r)
            
            fullmesh = UnitTriangleMesh(24,24);
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
        
        function result = solveForRadius(obj, r)
            % Resuelve el problema elástico para un radio dado usando Training
            % r: Radio de la inclusión
            % idx: Índice del radio (para logging)
            %
            % NOTA: Usa el mesh de referencia fijo para garantizar dimensiones consistentes
            % El material se calcula según el radio, pero el mesh permanece constante
                        
            meshRef = obj.createMesh(r);
            trainingData = Training(meshRef);  % Usa mesh de referencia fijo con material variable (r)
            
            % Obtain Kcoarse and Mcoasrse like for equilibrium problem
            processor = OfflineDataProcessor(trainingData);
            EIFEoper = processor.computeROMbasis(r);
            Kcoarse = EIFEoper.Kcoarse;
            Mcoarse = EIFEoper.Mcoarse;
            T = EIFEoper.T;
            
            
            result.r = r;
            result.u = T;
            result.Kcoarse = Kcoarse;
            result.Mcoarse = Mcoarse;
            result.mesh = trainingData.mesh;
        end
        
        function processTData(obj)
            % Procesa y formatea datos de T para CSV
            % Formato: [r, x, y, Tx1, Ty1, Tx2, Ty2, ..., Tx8, Ty8]
             nnodes = obj.allResults{1}.mesh.nnodes;
            TData = [];
            T_svd = zeros(nnodes*16,size(obj.radii,2));
            
            for j = 1:length(obj.allResults)
                result = obj.allResults{j};
                r = result.r;
                u = result.u; 
                mesh = result.mesh;
                
                % Reshape u: cada columna tiene [u_x1, u_y1, u_x2, u_y2, ...]
                % Necesitamos [Tx1, Ty1, Tx2, Ty2, ...] por nodo
                nnodes = mesh.nnodes;
                ndim = mesh.ndim;
                
                t_reshaped = zeros(nnodes, 16);  % 8 modos × 2 componentes
                for mode = 1:8
                    u_mode = u(:, mode);  % [nDofs × 1]
                    u_reshaped = reshape(u_mode, ndim, nnodes)';  % [nnodes × 2]
                    t_reshaped(:, 2*mode-1:2*mode) = u_reshaped;  % [Tx_mode, Ty_mode]
                end
                
                % COmbine to [r, x, y, Tx1, Ty1, ..., Tx8, Ty8]
                t_aux = [r * ones(nnodes, 1), mesh.coord, t_reshaped];
                TData = [TData; t_aux];
                %T_svd(:,j) = u(:);
            end
            
            obj.TData = TData; %used to train NN for Method A
            %obj.T_SVD = T_svd; %used for method B and C.
        end
        
        function processKMData(obj)
            % Procesa y formatea datos de K para CSV
            % format: [r, K11, K12, K13, ..., K88]
            
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
                       
           %
            fprintf('  Aplicando SVD...\n');
            [U_full, S_full, V_full] = svd(obj.T_SVD, 'econ'); %apply SVD
            
            tol = 1e-6;
            obj.k = sum(diag(S_full) > tol);  % Count significant singular values
           
                       
            % Truncar a k modos - for now we'll comment this and use the
            % full results
            obj.U = U_full; %(:, 1:obj.k);
            obj.S = S_full; %(1:obj.k, 1:obj.k);
            obj.V = V_full; %(:, 1:obj.k);
            fprintf('SVD completed. Modos retenidos: %d/%d\n', obj.k, min(size(obj.T_SVD)));
        end
        
        function processSVDData(obj)
            % Procesa datos SVD para entrenamiento de redes neuronales
            % Genera VTrainingData: [r, V1, V2, ..., Vk]
            % Each row is for a radius, each column for a mode

            % add radius column for training inoput
            obj.VTrainingData = [obj.radii(:), obj.V];
            
            fprintf('  Datos de entrenamiento V generados: [%d × %d]\n', ...
                size(obj.VTrainingData, 1), size(obj.VTrainingData, 2));
        end
                
    end
    
end
