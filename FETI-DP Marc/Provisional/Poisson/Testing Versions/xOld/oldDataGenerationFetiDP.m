classdef DataGenerationFetiDP < handle
    % DATAGENERATIONFETIDP Classe per resoldre un problema d'Elasticitat 2D
    % utilitzant FETI-DP i comparar-lo amb un mètode monolític.
    
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
        function obj = DataGenerationFetiDP()
            close all;
            obj.numSubdomains = [2 2]; 
            obj.nodeTol = 1e-10;
            obj.dofsPerNode = 2;
            
            % Creació de Malla Global i Subdominis
            referenceMesh = obj.createStructuredMesh();
            s.nsubdomains   = obj.numSubdomains;
            s.meshReference = referenceMesh;
            s.tolSameNode   = obj.nodeTol;
            m = MeshCreatorFromRVE.create(s);
            [obj.globalMesh, obj.localMeshes, ~, ~, ~, ~, ~] = m.create();
            
            disp(['Total elementos malla global: ', num2str(size(obj.globalMesh.connec, 1))]);
            disp(['Suma elementos subdominios: ', num2str(sum(cellfun(@(m) size(m.connec, 1), obj.localMeshes)))]);
            
            %obj.createBoundaryConditions(obj.globalMesh);
            
            % Identificar Dofs de Dirichlet (Només paret esquerra)
            obj.computeDirichletDofs();
            
            % Càlcul de Matrius Locals 
            obj.computeLocalMatrices();
            
            % Solver FETI-DP
            obj.fetiSolver = FetiDP(obj.globalMesh, obj.localMeshes, ...
                                    obj.localStiffness, obj.localForces, ...
                                    obj.nodeTol, obj.dofsPerNode, obj.dirichletDofs);
                                    
            obj.fetiSolver.visualizeFetiNodes();
            
            % Ensamblatge de FETI-DP
            [fMat, dBar] = obj.fetiSolver.assembleProblem();

            % Resolució a muerte sense precondicionador
            lambdaSol = fMat \ dBar; 
            
            % Resolució amb Solver Iteratiu PCG
            % tol = 1e-12;
            % maxit = 1000;
            % dirichletPrecond = @(r) obj.fetiSolver.applyDirichletPrecond(r);
            % 
            % tic
            % [lambdaSol, flag, ~, iter, resVec] = pcg(fMat, dBar, tol, maxit, dirichletPrecond);
            % toc
            % 
            % if flag == 0
            %     disp(['Convergència aconseguida en ', num2str(iter), ' iteracions!']);
            % end
            
            % Reconstrucció i Visualització de FETI-DP
            uFeti = obj.fetiSolver.reconstructGlobalSolution(lambdaSol, obj.globalMesh.nnodes);
            
            scaleFactor = 5; 
            obj.visualizeDeformedMesh(uFeti, scaleFactor, 'Solució FETI-DP');
            maxUyFeti = max(abs(uFeti(2:2:end)));
            fprintf('Màxim desplaçament en Y (FETI-DP) = %.4e\n', maxUyFeti);
            
            uMono = obj.solveMonolithicDirect();
            scaleFactor = 5; 
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
            % COMPUTEDIRICHLETDOFS Fixa la paret esquerra per a Elasticitat
            coords = obj.globalMesh.coord;
            minX = min(coords(:,1));
            
            isDirNode = abs(coords(:,1) - minX) < obj.nodeTol;
            dirNodes  = find(isDirNode);
            
            dDofs = [];
            for d = 1:obj.dofsPerNode
                dDofs = [dDofs; (dirNodes - 1) * obj.dofsPerNode + d];
            end
            obj.dirichletDofs = sort(dDofs);
        end
        
        function uDirect = solveMonolithicDirect(obj)
            % SOLVEMONOLITHICDIRECT Resol el problema global d'una sola vegada            
            uGlobal = LagrangianFunction.create(obj.globalMesh, obj.globalMesh.ndim, 'P1');
            material = obj.createMaterial(obj.globalMesh);
            
            cMat = material;
            weakForm = @(u,v) DDP(SymGrad(v), DDP(cMat, SymGrad(u)));
            kGlobal = IntegrateLHS(weakForm, uGlobal, uGlobal, obj.globalMesh, 'Domain', 2);
            
            totalDofs = obj.globalMesh.nnodes * obj.globalMesh.ndim;
            fGlobal   = zeros(totalDofs, 1);
            
            maxX = max(obj.globalMesh.coord(:,1));
            loadValue = -2; 
            
            for j = 1:obj.globalMesh.nnodes
                if abs(obj.globalMesh.coord(j,1) - maxX) < obj.nodeTol
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
            x1 = linspace(0, 1, 3);
            x2 = linspace(0, 1, 3);
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
            
            maxX = max(obj.globalMesh.coord(:,1));
            loadValue = -2; 
            
            for i = 1:numSub

                localMesh = obj.localMeshes{i};

                localMesh.plot();
                hold on;
                               
                uLoc = LagrangianFunction.create(localMesh, localMesh.ndim, 'P1');
                matLoc = obj.createMaterial(localMesh);
                
                cMat = matLoc;
                weakForm = @(u,v) DDP(SymGrad(v), DDP(cMat, SymGrad(u)));
                kCell{i} = IntegrateLHS(weakForm, uLoc, uLoc, localMesh, 'Domain', 2);
                
                fLocVec = zeros(uLoc.nDofs, 1);
            
                localCoords = localMesh.coord;
                for j = 1:size(localCoords, 1)
                    if abs(localCoords(j,1) - maxX) < obj.nodeTol
                        yDof = (j - 1) * localMesh.ndim + 2; 
                        fLocVec(yDof) = fLocVec(yDof) + loadValue;
                    end
                end
                
                fCell{i} = fLocVec;
            end
            
            obj.localStiffness = kCell;
            obj.localForces = fCell;
        end

        function [Dir,PL] = createRawBoundaryConditions(obj)
            minx = min(obj.globalMesh.coord(:,1));
            maxx = max(obj.globalMesh.coord(:,1));
            miny = min(obj.globalMesh.coord(:,2));
            maxy = max(obj.globalMesh.coord(:,2));
            tolBound = obj.nodeTol;
            isLeft   = @(coor) (abs(coor(:,1) - minx)   < tolBound);
            isRight  = @(coor) (abs(coor(:,1) - maxx)   < tolBound);
            isBottom = @(coor) (abs(coor(:,2) - miny)   < tolBound);
            isTop    = @(coor) (abs(coor(:,2) - maxy)   < tolBound);
            %             isMiddle = @(coor) (abs(coor(:,2) - max(coor(:,2)/2)) == 0);
            Dir{1}.domain    = @(coor) isLeft(coor);%| isRight(coor) ;
            Dir{1}.direction = [1,2];
            Dir{1}.value     = 0;

                        Dir{2}.domain    = @(coor) isRight(coor) ;
                        Dir{2}.direction = [1,2];
                        Dir{2}.value     = 0;

            PL.domain    = @(coor) isTop(coor);
            PL.direction = [2];
            PL.value     = [-0.1];
%                         PL.domain    = @(coor) isRight(coor);
%                         PL.direction = [1];
%                         PL.value     = [0.1];
        end

        function [bc,Dir,PL] = createBoundaryConditions(obj,mesh)
            [Dir,PL]  = obj.createRawBoundaryConditions();
            dirichletFun = [];
            for i = 1:numel(Dir)
                dir = DirichletCondition(obj.globalMesh, Dir{i});
                dirichletFun = [dirichletFun, dir];
            end

            pointload = PointLoad(mesh,PL);
            % need this because force applied in the face not in a point
            pointload.values        = pointload.values/size(pointload.dofs,1);
            fvalues                 = zeros(mesh.nnodes*mesh.ndim,1);
            fvalues(pointload.dofs) = pointload.values;
            fvalues                 = reshape(fvalues,mesh.ndim,[])';
            pointload.fun.setFValues(fvalues);

            s.pointloadFun = pointload;
            s.dirichletFun = dirichletFun;
            s.periodicFun  =[];
            s.mesh         = mesh;
            bc             = BoundaryConditions(s);
        end
    end
end