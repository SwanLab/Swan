classdef Tutorial_TopOpt2DLatticeXiMacroNN < handle

    properties (Access = private)
        mesh
        filters
        designVariable
        C
        dC
        physicalProblem
        compliance
        cost
        optimizer
    end

    methods (Access = public)

        function obj = Tutorial_TopOpt2DLatticeXiMacroNN()

            obj.init();
            obj.createMesh();
            obj.createDesignVariable();
            obj.createFilters();
            obj.createMaterial();
            obj.createElasticProblem();
            obj.createCompliance();
            obj.createCost();
            obj.createOptimizer();

        end

        function saveDesign(obj,fileName)

            fields = obj.designVariable.fields;
            [rho,~,mU] = obj.createDesignParameters();

            design.coord  = obj.mesh.coord;
            design.connec = obj.mesh.connec;
            design.rho    = rho*ones(obj.mesh.nnodes,1);
            design.a      = obj.filters{1}.compute(fields{1},2).fValues;
            design.xi     = obj.filters{2}.compute(fields{2},2).fValues;
            design.mU     = mU;

            save(fileName,'-struct','design');

        end

        function plotFields(obj)

            plotLatticeFields( ...
                obj.designVariable.fields, ...
                {'a','\xi'}, ...
                obj.createFieldLimits());

        end

    end

    methods (Access = private)

        function init(obj)

            close all;

        end

        function createMesh(obj)

            obj.mesh = TriangleMesh(2,1,100,50);

        end

        function createDesignVariable(obj)

            sA.fHandle = @(x) zeros(size(x(1,:,:)));
            sA.ndimf   = 1;
            sA.mesh    = obj.mesh;

            aAnalytical = AnalyticalFunction(sA);
            a = aAnalytical.project('P1');

            sXi.fHandle = @(x) pi/4*ones(size(x(1,:,:)));
            sXi.ndimf   = 1;
            sXi.mesh    = obj.mesh;

            xiAnalytical = AnalyticalFunction(sXi);
            xi = xiAnalytical.project('P1');

            sD.fields   = {a,xi};
            sD.type     = 'MultiField';
            sD.plotting = true;

            obj.designVariable = DesignVariable.create(sD);

        end

        function createFilters(obj)

            h = obj.mesh.computeMeanCellSize();

            sA.filterType   = 'PDE';
            sA.mesh         = obj.mesh;
            sA.boundaryType = 'Neumann';
            sA.metric       = 'Isotropy';
            sA.trial        = LagrangianFunction.create( ...
                obj.mesh,1,'P1');

            filterA = Filter.create(sA);
            filterA.updateEpsilon(6*h);

            sXi.filterType   = 'PDE';
            sXi.mesh         = obj.mesh;
            sXi.boundaryType = 'Neumann';
            sXi.metric       = 'Isotropy';
            sXi.trial        = LagrangianFunction.create( ...
                obj.mesh,1,'P1');

            filterXi = Filter.create(sXi);
            filterXi.updateEpsilon(6*h);

            obj.filters = {filterA,filterXi};

        end

        function createMaterial(obj)

            s.fileName      = 'latticeXiNN_rho08_Adam_64x3_cholesky_100k';
            s.mesh          = obj.mesh;
            s.young         = 1;
            s.parameterNames = {'a','xi'};

            hm = LatticeHomogenizedMaterialsReader(s);

            obj.C = @(designVariables) ...
                hm.obtainTensor(designVariables);

            obj.dC = @(designVariables) ...
                hm.obtainTensorDerivatives(designVariables);

        end

        function createElasticProblem(obj)

            s.mesh               = obj.mesh;
            s.scale              = 'MACRO';
            s.material           = [];
            s.dim                = '2D';
            s.boundaryConditions = obj.createBoundaryConditions();
            s.interpolationType  = 'LINEAR';
            s.solverType         = 'REDUCED';
            s.solverMode         = 'DISP';
            s.solverCase         = DirectSolver();

            obj.physicalProblem = ElasticProblem(s);

        end

        function c = createComplianceFromConstitutive(obj)

            s.mesh         = obj.mesh;
            s.stateProblem = obj.physicalProblem;

            c = ComplianceFromConstitutiveTensor(s);

        end

        function createCompliance(obj)

            s.mesh                       = obj.mesh;
            s.filters                    = obj.filters;
            s.complainceFromConstitutive = ...
                obj.createComplianceFromConstitutive();
            s.C                          = obj.C;
            s.dC                         = obj.dC;

            obj.compliance = ComplianceFunctional(s);

        end

        function createCost(obj)

            s.shapeFunctions{1} = obj.compliance;
            s.weights           = 1;
            s.Msmooth           = obj.createMassMatrix();

            obj.cost = Cost(s);

        end

        function M = createMassMatrix(obj)

            test = LagrangianFunction.create( ...
                obj.mesh,1,'P1');

            trial = LagrangianFunction.create( ...
                obj.mesh,1,'P1');

            MSingle = IntegrateLHS( ...
                @(u,v) DP(v,u), ...
                test, ...
                trial, ...
                obj.mesh, ...
                'Domain');

            MSingle = diag(sum(MSingle,1));

            nFields = obj.designVariable.obtainNumberVariables();

            M = kron(speye(nFields),MSingle);

        end

        function [rho,aMax,mU] = createDesignParameters(obj)

            rho  = 0.8;
            aMax = 0.8;
            mU   = latticeXiSettings().hole.mU;

        end

        function limits = createFieldLimits(obj)

            [rho,aMax,mU] = obj.createDesignParameters();

            [xiMin,xiMax] = XiLatticeDatasetGenerator.computeXiBounds( ...
                sqrt(1-rho),mU);

            limits = {[-aMax,aMax],[xiMin,xiMax]};

        end

        function [lb,ub] = createBounds(obj)

            fields = obj.designVariable.obtainDomainFunction();
            limits = obj.createFieldLimits();

            lb = [];
            ub = [];

            for iField = 1:numel(fields)

                n  = numel(fields{iField}.fValues);
                lb = [lb; limits{iField}(1)*ones(n,1)];
                ub = [ub; limits{iField}(2)*ones(n,1)];

            end

        end

        function createOptimizer(obj)

            [lb,ub] = obj.createBounds();

            s.monitoring      = true;
            s.cost            = obj.cost;
            s.designVariable  = obj.designVariable;
            s.maxIter         = 1000;
            s.lb              = lb;
            s.ub              = ub;

            opt = OptimizerProjectedGradient(s);

            opt.solveProblem();

            obj.optimizer = opt;

        end

        function bc = createBoundaryConditions(obj)

            xMax = max(obj.mesh.coord(:,1));
            yMax = max(obj.mesh.coord(:,2));

            isDir = @(coor) ...
                abs(coor(:,1)) == 0;

            isForce = @(coor) ...
                (abs(coor(:,1)) == xMax & ...
                abs(coor(:,2)) >= 0.4*yMax & ...
                abs(coor(:,2)) <= 0.6*yMax);

            sDir{1}.domain    = @(coor) isDir(coor);
            sDir{1}.direction = [1,2];
            sDir{1}.value     = 0;

            sPL{1}.domain    = @(coor) isForce(coor);
            sPL{1}.direction = 2;
            sPL{1}.value     = -1;

            dirichletFun = [];

            for i = 1:numel(sDir)

                dir = DirichletCondition( ...
                    obj.mesh, ...
                    sDir{i});

                dirichletFun = [dirichletFun,dir];

            end

            s.dirichletFun = dirichletFun;

            pointloadFun = [];

            for i = 1:numel(sPL)

                pl = TractionLoad( ...
                    obj.mesh, ...
                    sPL{i}, ...
                    'DIRAC');

                pointloadFun = [pointloadFun,pl];

            end

            s.pointloadFun = pointloadFun;
            s.periodicFun  = [];
            s.mesh         = obj.mesh;

            bc = BoundaryConditions(s);

        end

    end

end