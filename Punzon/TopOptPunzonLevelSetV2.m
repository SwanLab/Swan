classdef TopOptPunzonLevelSetV2 < handle

    properties (Access = private)
        mesh
        filter
        designVariable
        materialInterpolator
        physicalProblem
        compliance
        volume
        cost
        constraint
        primalUpdater
        optimizer
        filename
    end

    methods (Access = public)

        function obj = TopOptPunzonLevelSetV2()
            obj.init()
            obj.createMesh();
            obj.createDesignVariable();
            obj.createFilter();
            obj.createMaterialInterpolator();
            obj.createElasticProblem();
            obj.createComplianceFromConstiutive();
            obj.createCompliance();
            obj.createVolumeConstraint();
            obj.createCost();
            obj.createConstraint();
            obj.createPrimalUpdater();
            obj.createOptimizer();

            % Save monitoring and desginVariable fValues
            figure(2)
            saveas(gcf,'Monitoring_PunzonLevelSet.fig');
            obj.designVariable.fun.print('fValues_PunzonLevelSet');
        end

    end

    methods (Access = private)

        function createMesh(obj)
            file = 'punzon';
            obj.filename = file;
            a.fileName = file;
            s = FemDataContainer(a);
            obj.mesh = s.mesh;
        end

        function createDesignVariable(obj)
            % s.type = 'Full';
            % g      = GeometricalFunction(s);
            % lsFun  = g.computeLevelSetFunction(obj.mesh);
            % s.fun  = lsFun;

            s.fHandle = @(x) ones(size(x(1,:,:)));
            s.ndimf   = 1;
            s.mesh    = obj.mesh;
            aFun      = AnalyticalFunction(s);


            s.fun     = aFun.project('P1');
            s.mesh = obj.mesh;
            s.type = 'LevelSet';
            s.plotting = true;


            zMin     = min(obj.mesh.coord(:,3));
            isBottom = @(x) abs(x(:,3) - zMin) < 1e-6; % cara llisa
            guide1   = @(x) x(:,2)<= 20.179;
            guide2   = @(x) x(:,2)>= 76.729;

            s.isFixed  = obj.computeFixedVolumeDomain(@(x) guide1(x) | guide2(x) | isBottom(x), s.type);
            ls     = DesignVariable.create(s);
            obj.designVariable = ls;
        end

        function createFilter(obj)
            s.filterType = 'LUMP';
            s.mesh  = obj.mesh;
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            f = Filter.create(s);
            obj.filter = f;
        end

        function createMaterialInterpolator(obj)
            E0 = 1e-3;
            nu0 = 1/3;
            ndim = obj.mesh.ndim;
            matA.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E0,nu0);
            matA.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E0,nu0,ndim);


            E1 = 1;
            nu1 = 1/3;
            matB.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E1,nu1);
            matB.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E1,nu1,ndim);

            s.interpolation  = 'SIMPALL';
            s.dim            = '3D';
            s.matA = matA;
            s.matB = matB;

            m = MaterialInterpolator.create(s);
            obj.materialInterpolator = m;
        end

        function createElasticProblem(obj)
            s.mesh = obj.mesh;
            s.scale = 'MACRO';
            s.material = obj.createMaterial();
            s.dim = '3D';
            s.boundaryConditions = obj.createBoundaryConditions();
            s.interpolationType = 'LINEAR';
            s.solverType = 'REDUCED';
            s.solverMode = 'DISP';
            s.solverCase = DirectSolver();
            fem = ElasticProblem(s);
            obj.physicalProblem = fem;
        end

        function c = createComplianceFromConstiutive(obj)
            s.mesh         = obj.mesh;
            s.stateProblem = obj.physicalProblem;
            c = ComplianceFromConstitutiveTensor(s);
        end

        function createCompliance(obj)
            s.mesh                       = obj.mesh;
            s.filter                     = obj.filter;
            s.complainceFromConstitutive = obj.createComplianceFromConstiutive();
            s.material                   = obj.createMaterial();
            c = ComplianceFunctional(s);
            obj.compliance = c;
        end

        function uMesh = createBaseDomain(obj)
            sG.type          = 'Full';
            sG.length        = 1;
            sG.xCoorCenter   = 1.5;
            sG.yCoorCenter   = 0.5;
            g                = GeometricalFunction(sG);
            lsFun            = g.computeLevelSetFunction(obj.mesh);
            levelSet         = lsFun.fValues;
            s.backgroundMesh = obj.mesh;
            s.boundaryMesh   = obj.mesh.createBoundaryMesh();
            uMesh            = UnfittedMesh(s);
            uMesh.compute(levelSet);
        end

        function createVolumeConstraint(obj)
            s.mesh   = obj.mesh;
            s.filter = obj.filter;
            s.test = LagrangianFunction.create(obj.mesh,1,'P1');
            s.volumeTarget = 0.4;
            s.uMesh = obj.createBaseDomain();
            v = VolumeConstraint(s);
            obj.volume = v;
        end

        function createCost(obj)
            s.shapeFunctions{1} = obj.compliance;
            s.weights           = 1;
            s.Msmooth           = obj.createMassMatrix();
            obj.cost            = Cost(s);
        end

        function M = createMassMatrix(obj)
            n = obj.mesh.nnodes;
            h = obj.mesh.computeMinCellSize();
            M = h^2*sparse(1:n,1:n,ones(1,n),n,n);
        end

        function createConstraint(obj)
            s.shapeFunctions{1} = obj.volume;
            s.Msmooth           = obj.createMassMatrix();
            obj.constraint      = Constraint(s);
        end

        function createPrimalUpdater(obj)
            s.mesh = obj.mesh;
            obj.primalUpdater = SLERP(s);
        end

        function createOptimizer(obj)
            s.monitoring     = true;
            s.cost           = obj.cost;
            s.constraint     = obj.constraint;
            s.designVariable = obj.designVariable;
            s.maxIter        = 20;
            s.tolerance      = 1e-8;
            s.constraintCase = {'EQUALITY'};
            s.primalUpdater  = obj.primalUpdater;
            s.etaNorm        = 0.1;
            s.etaNormMin     = 0.05;
            s.gJFlowRatio    = 0.8;
            s.etaMax         = 10;
            s.etaMaxMin      = 0.5;
            %s.type           = '0';
            s.gif = false;
            s.gifName = [];
            s.printing = true;
            s.printName = 'PunzonLevelSet';
            opt = OptimizerNullSpace(s);
            opt.solveProblem();
            obj.optimizer = opt;
        end

        function m = createMaterial(obj)
            x = obj.designVariable;
            f = x.obtainDomainFunction();
            f = f{1}.project('P1');
            % f = obj.filter.compute(f{1},1);            
            s.type                 = 'DensityBased';
            s.density              = f;
            s.materialInterpolator = obj.materialInterpolator;
            s.dim                  = '3D';
            s.mesh                 = obj.mesh;
            m = Material.create(s);
        end

        function bc = createBoundaryConditions(obj)
            zMin = min(obj.mesh.coord(:,3));
            zMax = max(obj.mesh.coord(:,3));

            isTop = @(coor) abs(coor(:,3) - zMax) < 1e-6; % tornillos
            isForce = @(x) (x(3,:,:) - zMin) < 1e-2;

            sDir{1}.domain    = @(coor) isTop(coor);
            sDir{1}.direction = [1,2,3];
            sDir{1}.value     = 0;

            sPL{1}.domain    = isForce;
            [bMesh,~]=obj.mesh.createSingleBoundaryMesh();
            sPL{1}.fun       = ConstantFunction.create([0,0,1],bMesh);


            dirichletFun = [];
            for i = 1:numel(sDir)
                dir = DirichletCondition(obj.mesh, sDir{i});
                dirichletFun = [dirichletFun, dir];
            end
            s.dirichletFun = dirichletFun;

            pointloadFun = [];
            for i = 1:numel(sPL)
                pl = TractionLoad(obj.mesh, sPL{i}, 'FUNCTION');
                pointloadFun = [pointloadFun, pl];
            end
            s.pointloadFun = pointloadFun;

            s.periodicFun = [];
            s.mesh        = obj.mesh;

            bc = BoundaryConditions(s);
        end


        function isFixed = computeFixedVolumeDomain(obj,cond,type)
            coor  = obj.mesh.coord;
            nodes = find(cond(coor));
            isFixed.nodes = nodes;
            switch type
                case 'Density'
                    values = ones(size(nodes));
                case 'LevelSet'
                    values = -ones(size(nodes));
            end
            isFixed.values = values;
        end
    end

    methods(Access=private,Static)
        function init()
            close all;
        end
    end
end