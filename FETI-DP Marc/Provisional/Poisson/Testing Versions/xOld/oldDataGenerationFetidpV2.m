classdef DataGenerationFetidpV2 < handle
    % DATAGENERATIONFETIDP Classe per resoldre un problema d'Elasticitat 2D
    
    properties (Access = private)
        globalMesh       
        localMeshes      
        numSubdomains    
        nodeTol          
        dofsPerNode      
        
        dirichletDofs    
        localStiffness   
        localForces      
        fetiSolver       
    end
    
    methods (Access = public)
        function obj = DataGenerationFetidpV2()
            close all;
            obj.numSubdomains = [10 4]; 
            obj.nodeTol = 1e-13;
            obj.dofsPerNode = 2;
            
            referenceMesh = obj.createStructuredMesh();
            s.nsubdomains   = obj.numSubdomains;
            s.meshReference = referenceMesh;
            s.tolSameNode   = obj.nodeTol;
            m = MeshCreatorFromRVE.create(s);
            [obj.globalMesh, obj.localMeshes, ~, ~, ~, ~, ~] = m.create();
            
            obj.computeDirichletDofs();
            obj.computeLocalMatrices();
            
            obj.fetiSolver = FetidpV2(obj.globalMesh, obj.localMeshes, ...
                                    obj.localStiffness, obj.localForces, ...
                                    obj.nodeTol, obj.dofsPerNode, obj.dirichletDofs);
                                    
            obj.fetiSolver.visualizeFetiNodes();
            
            [fMat, dBar] = obj.fetiSolver.assembleProblem();
            
            warning('off', 'MATLAB:singularMatrix');
            warning('off', 'MATLAB:nearlySingularMatrix');
            lambdaSol = fMat \ dBar; 
            warning('on', 'MATLAB:singularMatrix');
            warning('on', 'MATLAB:nearlySingularMatrix');
            
            uFeti = obj.fetiSolver.reconstructGlobalSolution(lambdaSol, obj.globalMesh.nnodes);
            
            scaleFactor = 5; 
            obj.visualizeDeformedMesh(uFeti, scaleFactor, 'Solució FETI-DP');
            maxUyFeti = max(abs(uFeti(2:2:end)));
            fprintf('Màxim desplaçament en Y (FETI-DP) = %.4e\n', maxUyFeti);
            
            uMono = obj.solveMonolithicDirect();
            obj.visualizeDeformedMesh(uMono, scaleFactor, 'Solució Directe');
            relError = norm(uFeti - uMono) / norm(uMono);
            fprintf('Error relatiu entre FETI-DP i Monolític: %e\n', relError);
            
            if relError < 1e-10
                disp('FETI-DP coincideix amb el solver directe.');
            end           
        end
    end
    
    methods (Access = private)
        function computeDirichletDofs(obj)
            %% Viga empotrada
            % coords = obj.globalMesh.coord;
            % minX = min(coords(:,1));
            % 
            % isDirNode = abs(coords(:,1) - minX) < obj.nodeTol;
            % dirNodes  = find(isDirNode);
            % 
            % dDofs = [];
            % for d = 1:obj.dofsPerNode
            %     dDofs = [dDofs; (dirNodes - 1) * obj.dofsPerNode + d];
            % end
            % obj.dirichletDofs = sort(dDofs);

            %% Viga flexió
            coords = obj.globalMesh.coord;
            minX = min(coords(:,1));
            maxX = max(coords(:,1));
            midX = (minX + maxX) / 2; 
            
            % CANVI: Fixem només els nodes on la coordenada X és el centre
            isDirNode = abs(coords(:,1) - midX) < obj.nodeTol;
            dirNodes  = find(isDirNode);
            
            dDofs = [];
            for d = 1:obj.dofsPerNode
                dDofs = [dDofs; (dirNodes - 1) * obj.dofsPerNode + d];
            end
            obj.dirichletDofs = sort(dDofs);
        end
        
        function uDirect = solveMonolithicDirect(obj)
            uGlobal = LagrangianFunction.create(obj.globalMesh, obj.globalMesh.ndim, 'P1');
            material = obj.createMaterial(obj.globalMesh);
            
            cMat = material;
            weakForm = @(u,v) DDP(SymGrad(v), DDP(cMat, SymGrad(u)));
            kGlobal = IntegrateLHS(weakForm, uGlobal, uGlobal, obj.globalMesh, 'Domain', 2);
            
            totalDofs = obj.globalMesh.nnodes * obj.globalMesh.ndim;
            fGlobal   = zeros(totalDofs, 1);
            
            %% Viga empotrada
            % maxX = max(obj.globalMesh.coord(:,1));
            % loadValue = -10; 
            % 
            % for j = 1:obj.globalMesh.nnodes
            %     if abs(obj.globalMesh.coord(j,1) - maxX) < obj.nodeTol
            %         yDof = (j - 1) * obj.globalMesh.ndim + 2;
            %         fGlobal(yDof) = fGlobal(yDof) + loadValue;
            %     end
            % end

            %% Viga Flexió
            minX = min(obj.globalMesh.coord(:,1));
            maxX = max(obj.globalMesh.coord(:,1));
            loadValue = 10; % CANVI: +2 per estirar cap amunt
            
            for j = 1:obj.globalMesh.nnodes
                xCoord = obj.globalMesh.coord(j,1);
                % CANVI: Apliquem la força si és a la punta esquerra O a la dreta
                if abs(xCoord - minX) < obj.nodeTol || abs(xCoord - maxX) < obj.nodeTol
                    yDof = (j - 1) * obj.globalMesh.ndim + 2;
                    fGlobal(yDof) = fGlobal(yDof) + loadValue;
                end
            end
            
            freeDofs = setdiff(1:totalDofs, obj.dirichletDofs);
            
            kRed = kGlobal(freeDofs, freeDofs);
            fRed = fGlobal(freeDofs);
            
            uRed = kRed \ fRed;
            
            uDirect = zeros(totalDofs, 1);
            uDirect(freeDofs) = uRed;
            
            maxUyMono = max(abs(uDirect(2:2:end)));
            fprintf('Màxim desplaçament en Y (Directe) = %.4e\n', maxUyMono);
        end
        
        function visualizeDeformedMesh(obj, uGlobal, scaleFactor, titleStr)
            coords = obj.globalMesh.coord;
            connec = obj.globalMesh.connec;
            ndim   = obj.globalMesh.ndim;
            numNodes = size(coords, 1);
        
            uReshaped = reshape(uGlobal, ndim, numNodes)';
            defCoords = coords + scaleFactor * uReshaped;
            dispMag = sqrt(uReshaped(:,1).^2 + uReshaped(:,2).^2);
        
            figure('Name', ['Malla Deformada: ', titleStr], 'Color', 'w');
            hold on; axis equal;
        
            patch('Faces', connec, 'Vertices', coords, ...
                  'FaceColor', 'none', 'EdgeColor', [0.8 0.8 0.8], ...
                  'LineStyle', '--', 'DisplayName', 'Malla Original');
        
            patch('Faces', connec, 'Vertices', defCoords, ...
                  'FaceVertexCData', dispMag, 'FaceColor', 'interp', ...
                  'EdgeColor', '#333333', 'LineWidth', 0.5, 'DisplayName', 'Malla Deformada');
        
            colormap(jet);
            c = colorbar;
            c.Label.String = 'Magnitud del Desplaçament ||u||';
            title(sprintf('%s (Factor Escala: %g)', titleStr, scaleFactor));
            xlabel('X'); ylabel('Y');
        end
        
        function mS = createStructuredMesh(~)
            x1 = linspace(0, 2, 5);
            x2 = linspace(0, 1, 5);
            [xv, yv] = meshgrid(x1, x2);
            [F, V] = mesh2tri(xv, yv, zeros(size(xv)), 'x');
            
            s.coord  = V(:,1:2);
            s.connec = F;
            s.interpType = 'LINEAR';
            mS = Mesh.create(s);
        end
        
        function material = createMaterial(obj, mesh)
            [young, poisson] = obj.computeElasticProperties(mesh);
            s.type    = 'ISOTROPIC';
            s.ptype   = 'ELASTIC';
            s.ndim    = mesh.ndim;
            s.young   = young;
            s.poisson = poisson;
            material  = Material.create(s);
        end
        
        function [young, poisson] = computeElasticProperties(~, mesh)
            eMod = 100000;
            nu = 0.3;
            ePstr  = eMod / (1 - nu^2);
            nuPstr = nu / (1 - nu);
            
            young   = ConstantFunction.create(ePstr, mesh);
            poisson = ConstantFunction.create(nuPstr, mesh);
        end
        
        function computeLocalMatrices(obj)
            numSub = prod(obj.numSubdomains);
            kCell = cell(numSub, 1);
            fCell = cell(numSub, 1);
            
            %% Viga Empotrada
            % maxX = max(obj.globalMesh.coord(:,1));
            % loadValue = -10; 

            %% Viga Flexió
            minX = min(obj.globalMesh.coord(:,1));
            maxX = max(obj.globalMesh.coord(:,1));
            loadValue = 10;
            
            nodeMultiplicity = zeros(obj.globalMesh.nnodes, 1);
            for i = 1:numSub
                [~, globalNodes] = ismembertol(obj.localMeshes{i}.coord, obj.globalMesh.coord, obj.nodeTol, 'ByRows', true);
                nodeMultiplicity(globalNodes) = nodeMultiplicity(globalNodes) + 1;
            end
            
            for i = 1:numSub
                localMesh = obj.localMeshes{i};
                               
                uLoc = LagrangianFunction.create(localMesh, localMesh.ndim, 'P1');
                matLoc = obj.createMaterial(localMesh);
                
                cMat = matLoc;
                weakForm = @(u,v) DDP(SymGrad(v), DDP(cMat, SymGrad(u)));
                kCell{i} = IntegrateLHS(weakForm, uLoc, uLoc, localMesh, 'Domain', 2);
                
                fLocVec = zeros(uLoc.nDofs, 1);
                
                localCoords = localMesh.coord;
                [~, globalNodes] = ismembertol(localCoords, obj.globalMesh.coord, obj.nodeTol, 'ByRows', true);
                %% Viga Empotrada 
                % for j = 1:size(localCoords, 1)
                %     if abs(localCoords(j,1) - maxX) < obj.nodeTol
                %         yDof = (j - 1) * localMesh.ndim + 2; 
                %         gNode = globalNodes(j);
                %         fLocVec(yDof) = fLocVec(yDof) + (loadValue / nodeMultiplicity(gNode));
                %     end
                % end

                %% Viga Flexió
                for j = 1:size(localCoords, 1)
                    xCoord = localCoords(j,1);
                    % CANVI: Busquem puntes esquerra i dreta del subdomini local
                    if abs(xCoord - minX) < obj.nodeTol || abs(xCoord - maxX) < obj.nodeTol
                        yDof = (j - 1) * localMesh.ndim + 2; 
                        gNode = globalNodes(j);
                        fLocVec(yDof) = fLocVec(yDof) + (loadValue / nodeMultiplicity(gNode));
                    end
                end
                
                fCell{i} = fLocVec;
            end
            
            obj.localStiffness = kCell;
            obj.localForces = fCell;
        end
    end
end