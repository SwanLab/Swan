classdef TutorialFirst < handle

    properties (Access = private)
        mesh
        filterCompliance
        designVariable
        materialInterpolator
        physicalProblem
        compliance
        cost
        optimizer
        filterRegularization
    end

    methods (Access = public)

        function obj = TutorialFirst()
            obj.init()
            obj.createMesh();
            obj.createDesignVariable();
            obj.createFilterRegularization();
            obj.createElasticProblem();
            obj.createBaseDomain();
            obj.createComplianceFromConstiutive();
            obj.createCost();
            obj.createOptimizer();
        end
        function plotMicrostructuresOnMesh(obj, nGridX, nGridY, scaleFactor)
            if nargin < 2; nGridX     = 30;  end
            if nargin < 3; nGridY     = 15;  end
            if nargin < 4; scaleFactor = 0.8; end

            x = obj.designVariable;
            f = x.obtainDomainFunction();
            b_raw = f{1}.project('P1');
            b_filtered_fun = obj.filterRegularization.compute(b_raw, 1);
            b_vals = b_filtered_fun.fValues;

            coord = obj.mesh.coord;
            xMin  = min(coord(:,1));
            xMax  = max(coord(:,1));
            yMin  = min(coord(:,2));
            yMax  = max(coord(:,2));

            figure('Position', [50, 50, 1400, 700], 'Color', 'white');
            axis equal; axis off; hold on;
            xlim([-0.2, xMax + 0.4]);
            ylim([-0.15, yMax + 0.15]);

            x_pts    = linspace(0.03*xMax, 0.97*xMax, nGridX);
            y_pts    = linspace(0.03*yMax, 0.97*yMax, nGridY);
            baseSize = (xMax / nGridX) * scaleFactor;

            for ix = 1:nGridX
                for iy = 1:nGridY
                    xc = x_pts(ix);
                    yc = y_pts(iy);

                    dist = sqrt((coord(:,1)-xc).^2 + (coord(:,2)-yc).^2);
                    [~, idx] = min(dist);
                    b_local  = max(-0.6, min(0.6, b_vals(idx)));

                    a_local = exp(b_local^2);
                    d_local = (1 + b_local^2) / a_local;

                    v1 = baseSize * [a_local, b_local];
                    v2 = baseSize * [b_local, d_local];

                    centro = [xc, yc] - 0.5*(v1 + v2);

                    px = centro(1) + [0, v1(1), v1(1)+v2(1), v2(1), 0];
                    py = centro(2) + [0, v1(2), v1(2)+v2(2), v2(2), 0];

                    fill(px, py, [0,0,0], 'FaceAlpha', 0.75, ...
                        'EdgeColor', 'k', 'LineWidth', 0.4);

                    hf = 0.30;
                    local_pts = [0.5-hf/2, 0.5-hf/2;
                        0.5+hf/2, 0.5-hf/2;
                        0.5+hf/2, 0.5+hf/2;
                        0.5-hf/2, 0.5+hf/2];

                    hole_x = centro(1) + local_pts(:,1)*v1(1) + local_pts(:,2)*v2(1);
                    hole_y = centro(2) + local_pts(:,1)*v1(2) + local_pts(:,2)*v2(2);

                    fill(hole_x, hole_y, 'w', 'EdgeColor', 'none');
                end
            end

            % ── Engaste à ESQUERDA (x=0, face inteira) ──────────────────
            plot([xMin, xMin], [yMin, yMax], 'k-', 'LineWidth', 4);
            for y = linspace(yMin, yMax, 8)
                plot([xMin-0.12, xMin], [y-0.04, y], 'k-', 'LineWidth', 1.5);
            end

            % ── Força à DIREITA, centro vertical (35% a 65% da altura) ──
            % direction=2 (vertical), value=-1 (para baixo)
            yForce = yMax / 2;
            yF_min = 0.35 * yMax;
            yF_max = 0.65 * yMax;

            % Linha a vermelho mostrando a zona de aplicação
            plot([xMax, xMax], [yF_min, yF_max], 'r-', 'LineWidth', 3);

            % Seta para BAIXO no centro da zona
            quiver(xMax, yForce, 0, -0.18, 0, ...
                'r', 'LineWidth', 3, 'MaxHeadSize', 2);
            text(xMax+0.05, yForce - 0.22, 'F', 'FontSize', 14, ...
                'FontWeight', 'bold', 'Color', 'r', ...
                'HorizontalAlignment', 'center');

            title('Cantilever Beam — Microestruturas optimizadas', ...
                'FontSize', 13, 'FontWeight', 'bold');
        end

    end

    methods (Access = private)

        function init(obj)
            close all;
        end

        function createMesh(obj)
            
            obj.mesh = TriangleMesh(2,1,200,200);
        end

        function createDesignVariable(obj)
            s.fHandle = @(x) zeros(size(x(1,:,:)));
            s.ndimf   = 1;
            s.mesh    = obj.mesh;
            aFun      = AnalyticalFunction(s);

            sD.fun      = aFun.project('P1');
            sD.mesh     = obj.mesh;
            sD.type     = 'Density';
            sD.plotting = true;
            rho        = DesignVariable.create(sD);

            obj.designVariable = rho;
        end

       
        function createFilterRegularization(obj)
            eOverhmin = 8;
            s.filterType    = 'PDE';
            s.mesh          = obj.mesh;
            s.boundaryType  = 'Neumann';  
            s.metric        = 'Isotropy';
            s.trial         = LagrangianFunction.create(obj.mesh, 1, 'P1');
            f = Filter.create(s);
            
            epsilon = eOverhmin * obj.mesh.computeMeanCellSize();
            f.updateEpsilon(epsilon);
            
            obj.filterRegularization = f;
        end

     

        

        function createElasticProblem(obj)
            s.mesh = obj.mesh;
            s.scale = 'MACRO';
            s.material = obj.createMaterial();
            s.dim = '2D';
            s.boundaryConditions = obj.createBoundaryConditions();
            s.interpolationType = 'LINEAR';
            s.solverType = 'REDUCED';
            s.solverMode = 'DISP';
            s.solverCase = DirectSolver();
            fem = ElasticProblem(s);
            obj.physicalProblem = fem;
        end

        function uMesh = createBaseDomain(obj)
            levelSet         = -ones(obj.mesh.nnodes,1);
            s.backgroundMesh = obj.mesh;
            s.boundaryMesh   = obj.mesh.createBoundaryMesh();
            uMesh = UnfittedMesh(s);
            uMesh.compute(levelSet);
        end

        function c = createComplianceFromConstiutive(obj)
            s.mesh         = obj.mesh;
            s.stateProblem = obj.physicalProblem;
            c = ComplianceFromConstitutiveTensor(s);
        end
        function c = createCompliance(obj)
            s.mesh                        = obj.mesh;
            s.filter                      = obj.filterRegularization;
            s.complainceFromConstitutive  = obj.createComplianceFromConstiutive();
            s.material                    = obj.createMaterial();
            c = ComplianceFunctional(s);
            
        end


        
        function createCost(obj)
            s.shapeFunctions{1} = obj.createCompliance();
            s.weights           = 1;
            s.Msmooth           = obj.createMassMatrix();
            obj.cost            = Cost(s);
        end

        function M = createMassMatrix(obj)
            test  = LagrangianFunction.create(obj.mesh,1,'P1');
            trial = LagrangianFunction.create(obj.mesh,1,'P1'); 
            M = IntegrateLHS(@(u,v) DP(v,u),test,trial,obj.mesh,'Domain');
        end


        function createOptimizer(obj)
            s.cost           = obj.cost;
            s.designVariable = obj.designVariable;
            s.monitoring     = true;
            s.lb             = -0.6;
            s.ub             = 0.6;
            s.maxIter        = 200;
            opt              = OptimizerProjectedGradient(s);
            opt.solveProblem();
            obj.optimizer = opt;
            

            

            
        end

    
        function m = createMaterial(obj)
            x = obj.designVariable;
            
            
            s.density  = x;
            s.type     = 'HomogenizedMicrostructure';
            s.mesh     = obj.mesh;
            s.young    = 1.0;
            s.fileName = 'HomogenizationLattice4';
            m = MaterialFactory.create(s);
        end

        function bc = createBoundaryConditions(obj)
            xMax    = max(obj.mesh.coord(:,1));
            yMax    = max(obj.mesh.coord(:,2));
            isDir   = @(coor)  abs(coor(:,1))==0;
            isForce = @(coor)  (abs(coor(:,1))==xMax & abs(coor(:,2))>=0.35*yMax & abs(coor(:,2))<=0.65*yMax);

            sDir{1}.domain    = @(coor) isDir(coor);
            sDir{1}.direction = [1,2];
            sDir{1}.value     = 0;

            sPL{1}.domain    = @(coor) isForce(coor);
            sPL{1}.direction = 2;
            sPL{1}.value     = -1;

            dirichletFun = [];
            for i = 1:numel(sDir)
                dir = DirichletCondition(obj.mesh, sDir{i});
                dirichletFun = [dirichletFun, dir];
            end
            s.dirichletFun = dirichletFun;

            pointloadFun = [];
            for i = 1:numel(sPL)
                pl = TractionLoad(obj.mesh, sPL{i}, 'DIRAC');
                pointloadFun = [pointloadFun, pl];
            end
            s.pointloadFun = pointloadFun;

            s.periodicFun  = [];
            s.mesh = obj.mesh;
            bc = BoundaryConditions(s);
        end
    end
end