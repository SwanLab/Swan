classdef FetiDPPoissonScalability < handle    
    properties (Access = private)
        subdomainCases
        localNodesPerSide
        nodeTol
    end
    
    methods (Access = public)
        function obj = FetiDPPoissonScalability()
            close all;
            
            % --- CONFIGURACIÓN DEL ESTUDIO ---
            % Lista de configuraciones de subdominios a estudiar [Nx, Ny]
            obj.subdomainCases = {[2 2], [10 10], [20 20], [30 30], [40 40]};
            
            % Densidad de la malla local (H/h fijo = Constante)
            obj.localNodesPerSide = 8; 
            obj.nodeTol = 1e-10;
            
            % Ejecutar el bucle de escalabilidad
            obj.runScalabilityLoop();
        end
    end
    
    methods (Access = private)
        
        function runScalabilityLoop(obj)
            numCases = length(obj.subdomainCases);
            results  = zeros(numCases, 4); % [Subdominios Totales, DOFs Interfaz, Iteraciones, Kappa]
            
            % Creamos la malla local de referencia (RVE) UNA SOLA VEZ.
            % Así garantizamos que la relación H/h es estrictamente constante.
            refMesh = obj.createStructuredMesh(obj.localNodesPerSide);
            
            fprintf('\n=================================================================\n');
            fprintf('   ESTUDIO DE ESCALABILIDAD DÉBIL (FETI-DP + DIRICHLET) \n');
            fprintf('   Malla local fija: %dx%d nodos  (H/h constante)\n', obj.localNodesPerSide, obj.localNodesPerSide);
            fprintf('=================================================================\n');
            fprintf(' Subdominios | DOFs Interfaz | Iteraciones | Número de Condición (κ) \n');
            fprintf('-----------------------------------------------------------------\n');
            
            for i = 1:numCases
                subs = obj.subdomainCases{i};
                nSub = prod(subs);
                
                % 1. Generación de Malla
                s.nsubdomains   = subs;
                s.meshReference = refMesh;
                s.tolSameNode   = obj.nodeTol;
                m = MeshCreatorFromRVE.create(s);
                [globalMesh, localMeshes, ~, ~, ~, ~, ~] = m.create();
                
                % 2. Condiciones de contorn (Empotrado en bordes globales)
                dirDofs = obj.computeDirichletDofs(globalMesh);
                
                % 3. Càlcul de Matrius Locals (SENSE construir la K Global)
                [kLoc, fLoc] = obj.computeLocalMatrices(localMeshes);
                
                % 4. Solver FETI-DP
                feti = FetiDPPoissonV2(globalMesh, localMeshes, kLoc, fLoc, obj.nodeTol, dirDofs);
                [fMat, dBar] = feti.assembleProblem();
                feti.visualizeFetiNodes();
                
                % 5. PCG con Precondicionador de Dirichlet
                x0 = zeros(size(dBar));
                lambdaExact = fMat \ dBar; % Usado solo para cálculo de residuo exacto en PCG
                Pdir = @(r) feti.applyDirichletPrecond(r);
                
                [~, resHist, ~, ~] = PCG.solve(@(x) fMat*x, dBar, x0, Pdir, 1e-10, lambdaExact);
                
                % 6. Extracción de Mètriques
                iter = length(resHist);
                dofs = length(dBar);
                
                % Cálculo del número de condición Kappa
                M = feti.buildPrecondMatrix();
                kappa = condest(M * fMat);
                
                % Guardamos e imprimimos
                results(i, :) = [nSub, dofs, iter, kappa];
                fprintf(' %2dx%-2d = %-3d | %13d | %11d | %15.4f \n', subs(1), subs(2), nSub, dofs, iter, kappa);
            end
            fprintf('=================================================================\n\n');
            
            % Dibujar resultados
            obj.plotScalabilityResults(results);
        end
        
        function [kCell, fCell] = computeLocalMatrices(~, localMeshes)
            numSub = numel(localMeshes);
            kCell = cell(numSub, 1);
            fCell = cell(numSub, 1);
            
            for i = 1:numSub
                localMesh = localMeshes{i};
                uLoc = LagrangianFunction.create(localMesh, 1, 'P1');
                
                % LHS (Rigidez)
                weakFormK = @(u,v) DP(Grad(v), Grad(u));
                kCell{i} = IntegrateLHS(weakFormK, uLoc, uLoc, localMesh, 'Domain', 2);
                
                % RHS (Fuerzas)
                weakFormM = @(u,v) DP(v, u);
                mLocal = IntegrateLHS(weakFormM, uLoc, uLoc, localMesh, 'Domain', 2);
                fCell{i} = mLocal * ones(uLoc.nDofs, 1);
            end
        end
        
        function dirDofs = computeDirichletDofs(obj, globalMesh)
            coords = globalMesh.coord;
            minX = min(coords(:,1)); maxX = max(coords(:,1));
            minY = min(coords(:,2)); maxY = max(coords(:,2));
            
            isDirNode = (abs(coords(:,1) - minX) < obj.nodeTol) | ...
                        (abs(coords(:,1) - maxX) < obj.nodeTol) | ...
                        (abs(coords(:,2) - minY) < obj.nodeTol) | ...
                        (abs(coords(:,2) - maxY) < obj.nodeTol);
                        
            dirDofs = find(isDirNode);
        end
        
        function mS = createStructuredMesh(~, nPerSide)
            x1 = linspace(0, 1/5, nPerSide);
            x2 = linspace(0, 1/5, nPerSide);
            [xv, yv] = meshgrid(x1, x2);
            [F, V] = mesh2tri(xv, yv, zeros(size(xv)), 'x');
            s.coord  = V(:,1:2);
            s.connec = F;
            s.interpType = 'LINEAR';
            mS = Mesh.create(s);
        end
        
        function plotScalabilityResults(~, results)
            nSubs = results(:, 1);
            iters = results(:, 3);
            kappa = results(:, 4);
            
            figure('Name', 'Weak Scalability FETI-DP', 'Color', 'w', 'Position', [100, 100, 900, 400]);
            
            % Subplot 1: Número de Condición (Kappa)
            subplot(1,2,1);
            plot(nSubs, kappa, '-o', 'LineWidth', 2, 'MarkerSize', 6, 'Color', [0.85 0.33 0.10], 'MarkerFaceColor', [0.85 0.33 0.10]);
            title('Evolución del Número de Condición (\kappa)', 'FontSize', 12);
            xlabel('Número Total de Subdominios (N)', 'FontSize', 11);
            ylabel('\kappa (M^{-1} F)', 'FontSize', 11);
            grid on; ax = gca; ax.GridAlpha = 0.3;
            % Forzamos que el eje Y no empiece en 0 si no es necesario, para ver bien la línea plana
            y_margin = max(kappa)*0.1;
            ylim([min(kappa)-y_margin, max(kappa)+y_margin]);
            
            % Subplot 2: Iteraciones
            subplot(1,2,2);
            plot(nSubs, iters, '-s', 'LineWidth', 2, 'MarkerSize', 6, 'Color', [0.00 0.45 0.74], 'MarkerFaceColor', [0.00 0.45 0.74]);
            title('Evolución de las Iteraciones PCG', 'FontSize', 12);
            xlabel('Número Total de Subdominios (N)', 'FontSize', 11);
            ylabel('Iteraciones', 'FontSize', 11);
            grid on; ax = gca; ax.GridAlpha = 0.3;
            ylim([0, max(iters)+5]);
        end
    end
end