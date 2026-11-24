classdef Tutorial_TopOpt2DLatticeXiTGradientCheckStart < handle

    properties (Access = private)
        mesh
        filters
        designVariable
        C
        dC
        physicalProblem
        compliance
        volume
        cost
        constraint
    end

    methods (Access = public)

        function obj = Tutorial_TopOpt2DLatticeXiTGradientCheckStart()

            obj.init();
            obj.createMesh();
            obj.createDesignVariable();
            obj.createFilters();
            obj.createMaterial();
            obj.createElasticProblem();
            obj.createCompliance();
            obj.createVolumeConstraint();
            obj.createCost();
            obj.createConstraint();
            obj.checkGradients();

        end

    end

    methods (Access = private)

        function init(obj)

            close all;

        end

        function createMesh(obj)

            obj.mesh = TriangleMesh(2,1,40,20);

        end

        function createDesignVariable(obj)

            rho = obj.createField(@(x) 0.98*ones(size(x(1,:,:))));

            a   = obj.createField(@(x) zeros(size(x(1,:,:))));

            t   = obj.createField(@(x) 0.5*ones(size(x(1,:,:))));

            sD.fields   = {rho,a,t};
            sD.type     = 'MultiField';
            sD.plotting = false;

            obj.designVariable = DesignVariable.create(sD);

        end

        function f = createField(obj,fHandle)

            s.fHandle = fHandle;
            s.ndimf   = 1;
            s.mesh    = obj.mesh;

            fAnalytical = AnalyticalFunction(s);
            f = fAnalytical.project('P1');

        end

        function createFilters(obj)

            h = obj.mesh.computeMeanCellSize();

            sRho.filterType = 'LUMP';
            sRho.mesh       = obj.mesh;
            sRho.trial      = LagrangianFunction.create( ...
                obj.mesh,1,'P1');

            filterRho = Filter.create(sRho);

            sA.filterType   = 'PDE';
            sA.mesh         = obj.mesh;
            sA.boundaryType = 'Neumann';
            sA.metric       = 'Isotropy';
            sA.trial        = LagrangianFunction.create( ...
                obj.mesh,1,'P1');

            filterA = Filter.create(sA);
            filterA.updateEpsilon(6*h);

            sT.filterType   = 'PDE';
            sT.mesh         = obj.mesh;
            sT.boundaryType = 'Neumann';
            sT.metric       = 'Isotropy';
            sT.trial        = LagrangianFunction.create( ...
                obj.mesh,1,'P1');

            filterT = Filter.create(sT);
            filterT.updateEpsilon(6*h);

            obj.filters = {filterRho,filterA,filterT};

        end

        function createMaterial(obj)

            s.fileName      = 'latticeXiTNN_Adam_64x3_cholesky';
            s.mesh          = obj.mesh;
            s.young         = 1;
            s.parameterNames = {'rho','a','t'};

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

        function uMesh = createBaseDomain(obj)

            levelSet = -ones(obj.mesh.nnodes,1);

            s.backgroundMesh = obj.mesh;
            s.boundaryMesh   = obj.mesh.createBoundaryMesh();

            uMesh = UnfittedMesh(s);

            uMesh.compute(levelSet);

        end

        function createVolumeConstraint(obj)

            s.mesh             = obj.mesh;
            s.test             = LagrangianFunction.create( ...
                obj.mesh,1,'P1');
            s.volumeTarget     = 0.4;
            s.volumeFieldIndex = 1;
            s.uMesh            = obj.createBaseDomain();

            obj.volume = VolumeConstraint(s);

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

        function createConstraint(obj)

            s.shapeFunctions{1} = obj.volume;
            s.Msmooth           = obj.createMassMatrix();

            obj.constraint = Constraint(s);

        end

        function checkGradients(obj)

            x0         = obj.designVariable.fun.fValues;
            directions = obj.createDirections();
            names      = {'rho','a','t'};
            steps      = [1e-2 1e-3 1e-4 1e-5];

            obj.evaluateAt(x0);

            gJ = obj.cost.gradient;
            gV = obj.constraint.gradient;

            fprintf('\n%-5s %-6s %8s %14s %14s %10s\n', ...
                'func','field','step','finite diff','analytic','rel error');

            for iField = 1:numel(directions)

                d = directions{iField};

                for h = steps

                    [Jp,Vp] = obj.evaluateAt(x0 + h*d);
                    [Jm,Vm] = obj.evaluateAt(x0 - h*d);

                    obj.printRow('J',names{iField},h,(Jp-Jm)/(2*h),gJ'*d);
                    obj.printRow('V',names{iField},h,(Vp-Vm)/(2*h),gV'*d);

                end

            end

            obj.designVariable.update(x0);

        end

        function [J,V] = evaluateAt(obj,x)

            obj.designVariable.update(x);

            obj.cost.computeFunctionAndGradient(obj.designVariable);
            obj.constraint.computeFunctionAndGradient(obj.designVariable);

            J = obj.cost.value;
            V = obj.constraint.value;

        end

        function directions = createDirections(obj)

            dRho = obj.createField(@(x) ...
                0.5 + 0.5*cos(pi*x(1,:,:)/2).*sin(pi*x(2,:,:) + 0.5));

            dA = obj.createField(@(x) ...
                0.5*sin(pi*x(1,:,:)).*cos(pi*x(2,:,:)/2));

            dT = obj.createField(@(x) ...
                cos(pi*x(1,:,:) + 0.3).*cos(pi*x(2,:,:)));

            values = {dRho.fValues(:),dA.fValues(:),dT.fValues(:)};
            nTotal = sum(cellfun(@numel,values));

            directions = cell(numel(values),1);
            offset     = 0;

            for i = 1:numel(values)
                d = zeros(nTotal,1);
                d(offset + (1:numel(values{i}))) = values{i};
                directions{i} = d;
                offset = offset + numel(values{i});
            end

        end

        function printRow(~,func,field,h,fd,an)

            err = abs(fd-an)/max(abs(fd),eps);

            fprintf('%-5s %-6s %8.0e %14.6e %14.6e %10.2e\n', ...
                func,field,h,fd,an,err);

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