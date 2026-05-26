classdef TutorialFirstDensityOnly < handle

    properties (Access = private)
        mesh
        filterLUMP         % filtro LUMP para ρ
        physicalProblem
        materialMicro
        compliance
        volume
        cost
        constraint
        optimizer
    end
    properties (Access = public)
        designVariable
    end

    methods (Access = public)

        function hist = getCostHistory(obj)
            hist = obj.optimizer.costHistory;
        end

        function obj = TutorialFirstDensityOnly()
            obj.init()
            obj.createMesh();
            obj.createDesignVariable();
            obj.createFilterLUMP();
            obj.createMaterial();
            obj.createElasticProblem();
            obj.createBaseDomain();
            obj.createComplianceFromConstitutive();
            obj.createVolumeConstraint();
            obj.createCost();
            obj.createConstraint();
            obj.createOptimizer();
        end

               
        function plotMicrostructuresOnMesh(obj, nGridX, nGridY, scaleFactor)
            if nargin < 2; nGridX = 40; end
            if nargin < 3; nGridY = 20; end
            if nargin < 4; scaleFactor = 0.7; end

           
            x = obj.designVariable;
            f = x.obtainDomainFunction();
            
            
            rho_raw = f{1};
            rho_fun = obj.filterLUMP.compute(rho_raw, 1);
            rho_vals = rho_fun.fValues;
            
           
            b_vals = zeros(size(rho_vals));
           

            coord = obj.mesh.coord;

            xMax = max(coord(:,1));
            yMax = max(coord(:,2));
            yMin = min(coord(:,2));
            xMin = min(coord(:,1));

            figure('Position', [50, 50, 1400, 700], 'Color', 'white');
            axis equal; axis off; hold on;
            xlim([-0.2, xMax + 0.4]);
            ylim([-0.15, yMax + 0.15]);

            
            x_pts = linspace(0.03*xMax, 0.97*xMax, nGridX);
            y_pts = linspace(0.03*yMax, 0.97*yMax, nGridY);
            baseSize = (xMax / nGridX) * scaleFactor;

            for ix = 1:nGridX
                for iy = 1:nGridY
                    xc = x_pts(ix);
                    yc = y_pts(iy);

                    dist = sqrt((coord(:,1)-xc).^2 + (coord(:,2)-yc).^2);
                    [~, idx] = min(dist);
                    b_local = max(-0.6, min(0.6, b_vals(idx)));
                    rho_local = max(0, min(1, rho_vals(idx)));

                   
                    a_local = 1;
                    d_local = 1;

                    v1 = baseSize * [a_local, b_local];
                    v2 = baseSize * [b_local, d_local];

                    centro = [xc, yc] - 0.5*(v1 + v2);

                    px = centro(1) + [0, v1(1), v1(1)+v2(1), v2(1), 0];
                    py = centro(2) + [0, v1(2), v1(2)+v2(2), v2(2), 0];

                    gray = 1 - rho_local;
                    col = [gray, gray, gray];

                    fill(px, py, col, 'FaceAlpha', 0.85, ...
                        'EdgeColor', 'k', 'LineWidth', 0.5);

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

            
            plot([xMin, xMin], [yMin, yMax], 'k-', 'LineWidth', 4);
            for y = linspace(yMin, yMax, 8)
                plot([xMin - 0.1, xMin], [y, y], 'k-', 'LineWidth', 1.5);
            end

            xForce = xMax;
            yForce = (yMin + yMax) / 2;
            quiver(xForce + 0.05, yForce, 0, -0.2, 0, ...
                'r', 'LineWidth', 3, 'MaxHeadSize', 2);
            text(xForce + 0.15, yForce, 'F', 'FontSize', 14, ...
                'FontWeight', 'bold', 'Color', 'r');

            colormap(gray);
            caxis([0, 1]);
            colorbar('Location', 'eastoutside', 'FontSize', 10);
            ylabel(colorbar, 'Densidade ρ', 'FontSize', 10);

            title('Microestruturas - b=0 fixo (célula isotrópica), ρ variável', ...
                'FontSize', 12, 'FontWeight', 'bold');
        end
    end

    methods (Access = private)

        function init(obj)
            close all;
        end

        function createMesh(obj)
            obj.mesh = TriangleMesh(2, 1, 100, 60);
        end

        function createDesignVariable(obj)
            % Apenas ρ (b = 0 fixo)
            s.fHandle = @(x) ones(size(x(1,:,:)));
            s.ndimf = 1;
            s.mesh = obj.mesh;
            aFun = AnalyticalFunction(s);
            sD.fun = aFun.project('P1');
            sD.mesh = obj.mesh;
            sD.type = 'Density';
            sD.plotting = true;
            obj.designVariable = DesignVariable.create(sD);
        end

        function createFilterLUMP(obj)
            s.filterType = 'LUMP';
            s.mesh = obj.mesh;
            s.trial = LagrangianFunction.create(obj.mesh, 1, 'P1');
            obj.filterLUMP = Filter.create(s);
        end

        function createMaterial(obj)
            % Material com b = 0 fixo
            s.type = 'HomogenizedMicrostructureBZero';
            s.mesh = obj.mesh;
            s.young = 1.0;
            s.fileName = 'HomogenizationLattice4';
            s.density = obj.designVariable;
            obj.materialMicro = MaterialFactory.create(s);
        end

        function createElasticProblem(obj)
            s.mesh = obj.mesh;
            s.scale = 'MACRO';
            s.material = obj.materialMicro;
            s.dim = '2D';
            s.boundaryConditions = obj.createBoundaryConditions();
            s.interpolationType = 'LINEAR';
            s.solverType = 'REDUCED';
            s.solverMode = 'DISP';
            s.solverCase = DirectSolver();
            obj.physicalProblem = ElasticProblem(s);
        end

        function uMesh = createBaseDomain(obj)
            levelSet = -ones(obj.mesh.nnodes, 1);
            s.backgroundMesh = obj.mesh;
            s.boundaryMesh = obj.mesh.createBoundaryMesh();
            uMesh = UnfittedMesh(s);
            uMesh.compute(levelSet);
        end

        function createComplianceFromConstitutive(obj)
            s.mesh = obj.mesh;
            s.stateProblem = obj.physicalProblem;
            obj.compliance = ComplianceFromConstitutiveTensor(s);
        end

        function c = createCompliance(obj)
            s.mesh = obj.mesh;
            s.filter = obj.filterLUMP;
            s.complainceFromConstitutive = obj.compliance;
            s.material = obj.materialMicro;
            c = ComplianceFunctional(s);
        end

        function createVolumeConstraint(obj)
            s.mesh = obj.mesh;
            s.filter = obj.filterLUMP;
            s.test = LagrangianFunction.create(obj.mesh, 1, 'P1');
            s.volumeTarget = 0.4;
            s.uMesh = obj.createBaseDomain();
            v = VolumeConstraint(s);
            obj.volume = v;
        end

        function createCost(obj)
            s.shapeFunctions{1} = obj.createCompliance();
            s.weights = 1;
            s.Msmooth = obj.createMassMatrix();
            obj.cost = Cost(s);
        end

        function createConstraint(obj)
            s.shapeFunctions{1} = obj.volume;
            s.Msmooth = obj.createMassMatrix();
            obj.constraint = Constraint(s);
        end

        function M = createMassMatrix(obj)
            test = LagrangianFunction.create(obj.mesh, 1, 'P1');
            trial = LagrangianFunction.create(obj.mesh, 1, 'P1');
            M = IntegrateLHS(@(u,v) DP(v,u), test, trial, obj.mesh, 'Domain');
        end

        function p = createPrimalUpdater(obj)
            n = obj.mesh.nnodes;
            s.lb = zeros(n, 1);
            s.ub = ones(n, 1);
            s.tauMax = 200;
            s.tau = [];
            p = ProjectedGradient(s);
        end

        function createOptimizer(obj)
            s.monitoring = true;
            s.cost = obj.cost;
            s.constraint = obj.constraint;
            s.designVariable = obj.designVariable;
            s.maxIter = 1000;
            s.tolerance = 1e-8;
            s.constraintCase = {'EQUALITY'};
            s.etaNorm = 0.01;
            s.gJFlowRatio = 2;
            s.primalUpdater = obj.createPrimalUpdater();
            s.gif = false;
            s.gifName = [];
            s.printing = false;
            s.printName = [];
            opt = OptimizerNullSpace(s);
            opt.solveProblem();
            obj.optimizer = opt;
        end

        function bc = createBoundaryConditions(obj)
            xMax = max(obj.mesh.coord(:,1));
            yMax = max(obj.mesh.coord(:,2));
            yMin = min(obj.mesh.coord(:,2));
            xMin = min(obj.mesh.coord(:,1));

            isDir = @(coor) abs(coor(:,1) - xMin) < 1e-12;
            isForce = @(coor) abs(coor(:,1) - xMax) < 1e-12 & ...
                              coor(:,2) >= 0.4*yMax & coor(:,2) <= 0.6*yMax;

            sDir{1}.domain = @(coor) isDir(coor);
            sDir{1}.direction = [1,2];
            sDir{1}.value = 0;

            sPL{1}.domain = @(coor) isForce(coor);
            sPL{1}.direction = 2;
            sPL{1}.value = -1;

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

            s.periodicFun = [];
            s.mesh = obj.mesh;
            bc = BoundaryConditions(s);
        end
    end
end