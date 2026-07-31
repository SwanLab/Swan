classdef TutorialFirstBOnly < handle

    properties (Access = public)
        designVariable
        bSmooth
    end

    properties (Access = private)
        mesh
        physicalProblem
        materialMicro
        compliance
        cost
        optimizer
        filterRegularization
        baseDomain
    end

    methods (Access = public)

        function obj = TutorialFirstBOnly()

            obj.init();
            obj.createMesh();
            obj.createDesignVariable();
            obj.createFilterRegularization();
            obj.createBaseDomain();
            obj.createMaterial();
            obj.createElasticProblem();
            obj.createComplianceFromConstitutive();
            obj.createCost();
            obj.createOptimizer();

            obj.bSmooth = obj.filterRegularization.compute( ...
                obj.designVariable.fun,3);

        end

        function hist = getCostHistory(obj)
            hist = obj.optimizer.costHistory;
        end

        function hist = getComplianceHistory(obj)
            hist = obj.optimizer.complianceHistory;
        end
        function beta = getBValues(obj)
            % Variavel de projeto final, antes do filtro
            beta = obj.designVariable.fun.fValues;
        end

        function betaSmooth = getBSmoothValues(obj)
            % Variavel de projeto final filtrada
            betaSmooth = obj.bSmooth.fValues;
        end

        function bGeomSmooth = getBGeomSmoothValues(obj)
            % Parametro geometrico fisico da microestrutura
            %
            % Parametrizacao usada:
            % bGeom = bMax*(1-(beta/bMax)^2)

            bMax = 0.6;

            betaSmooth = obj.bSmooth.fValues;

            bGeomSmooth = bMax.*( ...
                1-(betaSmooth./bMax).^2);

            % Protecao contra pequenos erros numericos
            bGeomSmooth = min( ...
                max(bGeomSmooth,0), ...
                bMax);
        end

        function meshData = getMeshData(obj)
            % Guarda somente os dados numericos necessarios da malha

            meshData.coord  = obj.mesh.coord;
            meshData.connec = obj.mesh.connec;
            meshData.ndim   = obj.mesh.ndim;
            meshData.nnodes = obj.mesh.nnodes;
        end

       

    end

    methods (Access = private)

        function init(obj)
            close all;
        end

        function createMesh(obj)

            
            obj.mesh = TriangleMesh(2,1,200,180);

        end

        function createDesignVariable(obj)

           
            s.fHandle = @(x) zeros(size(x(1,:,:)));
            s.ndimf   = 1;
            s.mesh    = obj.mesh;

            aFun = AnalyticalFunction(s);
            funB = aFun.project('P1');

            sD.fun      = funB;
            sD.mesh     = obj.mesh;
            sD.type     = 'Density';
            sD.plotting = true;

            obj.designVariable = DesignVariable.create(sD);

        end

        function createFilterRegularization(obj)

            eOverhmin = 6;

            s.filterType   = 'PDE';
            s.mesh         = obj.mesh;
            s.boundaryType = 'Neumann';
            s.metric       = 'Isotropy';
            s.trial        = LagrangianFunction.create( ...
                obj.mesh,1,'P1');

            f = Filter.create(s);

            epsilon = eOverhmin*obj.mesh.computeMeanCellSize();
            f.updateEpsilon(epsilon);

            obj.filterRegularization = f;

        end

        function createBaseDomain(obj)

            levelSet = -ones(obj.mesh.nnodes,1);

            s.backgroundMesh = obj.mesh;
            s.boundaryMesh   = obj.mesh.createBoundaryMesh();

            uMesh = UnfittedMesh(s);
            uMesh.compute(levelSet);

            obj.baseDomain = uMesh;

        end

        function createMaterial(obj)

            
            s.type     = 'HomogenizedMicrostructureBOnly';
            s.mesh     = obj.mesh;
            s.young    = 1.0;
            s.fileName = 'Homogenization_bOnly_rho040';
            s.density  = obj.designVariable;

            obj.materialMicro = MaterialFactory.create(s);

        end

        function createElasticProblem(obj)

            s.mesh               = obj.mesh;
            s.scale              = 'MACRO';
            s.material           = obj.materialMicro;
            s.dim                = '2D';
            s.boundaryConditions = obj.createBoundaryConditions();
            s.interpolationType  = 'LINEAR';
            s.solverType         = 'REDUCED';
            s.solverMode         = 'DISP';
            s.solverCase         = DirectSolver();

            obj.physicalProblem = ElasticProblem(s);

        end

        function createComplianceFromConstitutive(obj)

            s.mesh         = obj.mesh;
            s.stateProblem = obj.physicalProblem;

            obj.compliance = ...
                ComplianceFromConstitutiveTensor(s);

        end

        function c = createCompliance(obj)

            s.mesh                       = obj.mesh;
            s.filter                     = obj.filterRegularization;
            s.complainceFromConstitutive = obj.compliance;
            s.material                   = obj.materialMicro;

            c = ComplianceFunctional(s);

        end

        function createCost(obj)

            s.shapeFunctions{1} = obj.createCompliance();
            s.weights           = 1;
            s.Msmooth           = obj.createMassMatrix();

            obj.cost = Cost(s);

        end

        function M = createMassMatrix(obj)

            test = LagrangianFunction.create( ...
                obj.mesh,1,'P1');

            trial = LagrangianFunction.create( ...
                obj.mesh,1,'P1');

            M = IntegrateLHS( ...
                @(u,v) DP(v,u), ...
                test,trial,obj.mesh,'Domain');

        end

        function createOptimizer(obj)

            s.cost           = obj.cost;
            s.designVariable = obj.designVariable;
            s.monitoring     = true;

            
            s.lb = -0.6;
            s.ub =  0.6;

            s.maxIter   = 1000;
            s.tolerance = 1e-8;
            
            

            opt = OptimizerProjectedGradient(s);
            opt.solveProblem();

            obj.optimizer = opt;

        end

        function bc = createBoundaryConditions(obj)

            xMin = min(obj.mesh.coord(:,1));
            xMax = max(obj.mesh.coord(:,1));
            yMax = max(obj.mesh.coord(:,2));

            isDir = @(coor) ...
                abs(coor(:,1)-xMin) < 1e-12;

            isForce = @(coor) ...
                abs(coor(:,1)-xMax) < 1e-12 & ...
                coor(:,2) >= 0.4*yMax & ...
                coor(:,2) <= 0.6*yMax;

            sDir{1}.domain    = @(coor) isDir(coor);
            sDir{1}.direction = [1,2];
            sDir{1}.value     = 0;

            sPL{1}.domain    = @(coor) isForce(coor);
            sPL{1}.direction = 2;
            sPL{1}.value     = -1;

            dirichletFun = [];

            for i = 1:numel(sDir)

                dirichletFun = [ ...
                    dirichletFun, ...
                    DirichletCondition(obj.mesh,sDir{i})];

            end

            pointloadFun = [];

            for i = 1:numel(sPL)

                pointloadFun = [ ...
                    pointloadFun, ...
                    TractionLoad( ...
                    obj.mesh,sPL{i},'DIRAC')];

            end

            s.dirichletFun = dirichletFun;
            s.pointloadFun = pointloadFun;
            s.periodicFun  = [];
            s.mesh         = obj.mesh;

            bc = BoundaryConditions(s);

        end

    end

end