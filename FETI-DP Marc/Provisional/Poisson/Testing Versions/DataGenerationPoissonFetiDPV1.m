classdef DataGenerationPoissonFetiDPV1 < handle
    % DATAGENERATIONPOISSONFETIDP Classe per resoldre l'equació de Poisson
    % utilitzant FETI-DP i comparar amb el mètode directe.
    
    properties (Access = private)
        globalMesh       
        localMeshes      
        numSubdomains    
        nodeTol          
        
        localStiffness   
        localForces      
        fetiSolver       
    end
    
    methods (Access = public)
        function obj = DataGenerationPoissonFetiDPV1()
            close all;
            obj.numSubdomains = [25 7]; 
            obj.nodeTol = 1e-10;
            
            % Creació de Malla Global 
            referenceMesh = obj.createStructuredMesh();
            s.nsubdomains   = obj.numSubdomains;
            s.meshReference = referenceMesh;
            s.tolSameNode   = obj.nodeTol;
            m = MeshCreatorFromRVE.create(s);
            [obj.globalMesh, obj.localMeshes, ~, ~, ~, ~, ~] = m.create();
            
            % Càlcul de Matrius Locals de Poisson
            obj.computeLocalPoissonMatrices();
            
            % Solver FETI-DP per Poisson
            obj.fetiSolver = FetiDPPoisson(obj.globalMesh, obj.localMeshes, ...
                                           obj.localStiffness, obj.localForces, ...
                                           obj.nodeTol);
                                    
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
            
            obj.visualizePoissonSolution(uFeti, 'Solució FETI-DP (Poisson)');
            maxUFeti = max(abs(uFeti));
            fprintf('Màxim valor u (FETI-DP) = %.4e\n', maxUFeti);
            
            % Comparació amb el Solver Monolític
            uMono = obj.solveMonolithicPoisson();

            obj.visualizePoissonSolution(uMono, 'Solució Monolítica (Poisson)');
            
            relError = norm(uFeti - uMono) / norm(uMono);
            fprintf('Error relatiu entre FETI-DP i Monolític: %e\n', relError);
            
            if relError < 1e-10
                disp('La solució FETI-DP coincideix amb el solver directe.');
            end           
        end
    end
    
    methods (Access = private)
        function uDirect = solveMonolithicPoisson(obj)
            % SOLVEMONOLITHICPOISSON Resol el problema de Poisson globalment
            disp('--- Ensamblant el Sistema Monolític Global (Poisson) ---');
            
            uGlobal = LagrangianFunction.create(obj.globalMesh, 1, 'P1');
            
            % Matriu de Rigidesa (Laplacià)
            weakFormK = @(u,v) DP(Grad(v), Grad(u));
            kGlobal = IntegrateLHS(weakFormK, uGlobal, uGlobal, obj.globalMesh, 'Domain', 2);
            
            % Càrrega volumètrica f=1 a tot el domini integrada amb la Matriu de Massa
            weakFormM = @(u,v) DP(v, u);
            mGlobal = IntegrateLHS(weakFormM, uGlobal, uGlobal, obj.globalMesh, 'Domain', 2);
            fGlobal = mGlobal * ones(uGlobal.nDofs, 1);
            
            % Identificar Dirichlet en tot el perímetre
            coords = obj.globalMesh.coord;
            minX = min(coords(:,1)); maxX = max(coords(:,1));
            minY = min(coords(:,2)); maxY = max(coords(:,2));
            
            isDirNode = (abs(coords(:,1) - minX) < obj.nodeTol) | ...
                        (abs(coords(:,1) - maxX) < obj.nodeTol) | ...
                        (abs(coords(:,2) - minY) < obj.nodeTol) | ...
                        (abs(coords(:,2) - maxY) < obj.nodeTol);
            dirichletDofs = find(isDirNode);
            
            freeDofs = setdiff(1:uGlobal.nDofs, dirichletDofs);
            
            disp('Resolent U = K \ F (Directe) ...');
            kRed = kGlobal(freeDofs, freeDofs);
            fRed = fGlobal(freeDofs);
            
            uRed = kRed \ fRed;
            
            uDirect = zeros(uGlobal.nDofs, 1);
            uDirect(freeDofs) = uRed;
            
            maxUyMono = max(abs(uDirect));
            fprintf('Màxim valor escalar (Monolític) = %.4e\n', maxUyMono);
        end
        
        function visualizePoissonSolution(obj, uGlobal, titleStr)
            % VISUALIZEPOISSONSOLUTION Dibuixa un mapa de calor del camp escalar
            coords = obj.globalMesh.coord;
            connec = obj.globalMesh.connec;
        
            figure('Name', titleStr, 'Color', 'w');
            hold on; axis equal;
        
            patch('Faces', connec, 'Vertices', coords, ...
                  'FaceVertexCData', uGlobal, 'FaceColor', 'interp', ...
                  'EdgeColor', '#333333', 'LineWidth', 0.5);
        
            colormap(jet);
            c = colorbar;
            c.Label.String = 'Camp Escalar u';
            title(titleStr);
            xlabel('X'); ylabel('Y');
        end
        
        function mS = createStructuredMesh(~)
            % CREATESTRUCTUREDMESH Crea la malla base
            x1 = linspace(0, 1/3, 3);
            x2 = linspace(0, 1/3, 3);
            [xv, yv] = meshgrid(x1, x2);
            [F, V] = mesh2tri(xv, yv, zeros(size(xv)), 'x');
            
            s.coord  = V(:,1:2);
            s.connec = F;
            s.interpType = 'LINEAR';
            mS = Mesh.create(s);
        end
        
        function computeLocalPoissonMatrices(obj)
            % COMPUTELOCALPOISSONMATRICES Ensambla el Laplacià i càrrega volumètrica f=1
            numSub = prod(obj.numSubdomains);
            kCell = cell(numSub, 1);
            fCell = cell(numSub, 1);
            
            for i = 1:numSub
                localMesh = obj.localMeshes{i};
                               
                uLoc = LagrangianFunction.create(localMesh, 1, 'P1');
                
                % Laplacià
                weakFormK = @(u,v) DP(Grad(v), Grad(u));
                kCell{i} = IntegrateLHS(weakFormK, uLoc, uLoc, localMesh, 'Domain', 2);
                
                % Matriu de massa per a la càrrega constant f=1.0
                weakFormM = @(u,v) DP(v, u);
                mLocal = IntegrateLHS(weakFormM, uLoc, uLoc, localMesh, 'Domain', 2);
                
                fCell{i} = mLocal * ones(uLoc.nDofs, 1);
            end
            
            obj.localStiffness = kCell;
            obj.localForces = fCell;
        end
    end
end