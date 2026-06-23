classdef LevelSetGripper6x6Thick < handle

    properties (Access = private)
        mesh
        filename
        filter
        designVariable
        materialInterpolator
        physicalProblem
        compliance
        penalty
        perimeter
        volume
        cost
        constraint
        primalUpdater
        optimizer
    end

    methods (Access = public)

        function obj = LevelSetGripper6x6Thick(hTar)
            obj.init()
            obj.createMesh();
            obj.createDesignVariable();
            obj.createFilter();
            obj.createMaterialInterpolator();
            obj.createElasticProblem();
            obj.createCompliance();
            obj.createPerimeterPenalty();
            obj.createThicknessConstraint(hTar);
            obj.createVolumeConstraint();
            obj.createCost();
            obj.createConstraint();
            obj.createPrimalUpdater();
            obj.createOptimizer();

            saveas(gcf,['Paper/Domain6x6/MonitoringLevelSetGripper6x6Thick',num2str(hTar),'.fig']);
            obj.designVariable.fun.print(['Paper/Domain6x6/LevelSetGripper6x6Thick',num2str(hTar),'fValues']);
        end

    end

    methods (Access = private)

        function init(obj)
            close all;
        end

        function createMesh(obj)
            file = 'Gripping';
            obj.filename = file;
            a.fileName = file;
            s = FemDataContainer(a);
            obj.mesh = s.mesh;
        end

        function createDesignVariable(obj)
            s.type             = 'Holes';
            s.dim              = 2;
            s.nHoles           = [50,50]; % 100, 100
            s.totalLengths     = [1,1];
            s.phiZero          = 0.4;
            s.phases           = [pi/2,0];
            g                  = GeometricalFunction(s);
            lsFun              = g.computeLevelSetFunction(obj.mesh);
            lsFun.setFValues(lsFun.fValues - 1);
            s.fun              = lsFun;
            s.mesh             = obj.mesh;
            s.type             = 'LevelSet';
            s.plotting         = false;
            ls                 = DesignVariable.create(s);
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
            s.dim            = '2D';
            s.matA = matA;
            s.matB = matB;

            m = MaterialInterpolator.create(s);
            obj.materialInterpolator = m;
        end

        function m = createMaterial(obj)
            f = obj.designVariable.fun;           
            s.type                 = 'DensityBased';
            s.density              = f;
            s.materialInterpolator = obj.materialInterpolator;
            s.dim                  = '2D';
            s.mesh                 = obj.mesh;
            m = Material.create(s);
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

        function createCompliance(obj)
            s.mesh         = obj.mesh;
            s.filter       = obj.filter;
            s.material     = obj.createMaterial();
            s.stateProblem = obj.physicalProblem;
            s.filename     = obj.filename;
            c = NonSelfAdjointComplianceFunctional(s);
            obj.compliance = c;
        end

        function createPerimeterPenalty(obj)
            sF.mesh       = obj.mesh;
            sF.filterType = 'PDE';
            sF.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
            f             = Filter.create(sF);

            h             = obj.mesh.computeMeanCellSize();
            s.mesh        = obj.mesh;
            s.uMesh       = obj.createBaseDomain();
            s.filter      = f;
            s.epsilon     = 3*h;
            s.value0      = 4;
            s.signInitial = -0.25;
            s.signFinal   = 0.05;
            s.tarVolume   = 0.5;
            obj.penalty   = InterfaceFunctional(s);
        end

        function uMesh = createBaseDomain(obj)
            levelSet         = -ones(obj.mesh.nnodes,1);
            s.backgroundMesh = obj.mesh;
            s.boundaryMesh   = obj.mesh.createBoundaryMesh();
            uMesh = UnfittedMesh(s);
            uMesh.compute(levelSet);
        end

        function uMesh = createBaseDomainIP(obj,x0,y0)
            s.type             = 'Rectangle';
            s.xCoorCenter      = x0;
            s.yCoorCenter      = y0;
            s.xSide            = 1/6;
            s.ySide            = 1/6;
            g                  = GeometricalFunction(s);
            lsFun              = g.computeLevelSetFunction(obj.mesh);
            sUm.backgroundMesh = obj.mesh;
            sUm.boundaryMesh   = obj.mesh.createBoundaryMesh();
            uMesh              = UnfittedMesh(sUm);
            uMesh.compute(lsFun.fValues);
        end

        function createThicknessConstraint(obj,hTar)
            sF.mesh        = obj.mesh;
            sF.filterType  = 'PDE';
            sF.boundaryType = 'Robin';
            sF.trial       = LagrangianFunction.create(obj.mesh,1,'P1');
            f              = Filter.create(sF);

            h         = obj.mesh.computeMeanCellSize();
            s.mesh    = obj.mesh;
            s.filter  = f;
            s.epsilon = 3*h;
            s.value0 = 1;
            s.tarVolume = 0.5;
            s.test = LagrangianFunction.create(obj.mesh,1,'P1');
            s.tau = 0.1/6;

            s.target = hTar;
            s.target0 = s.target/100;

            x0 = repmat([1/12,1/4,5/12,7/12,3/4,11/12],[1,6]);
            y0 = [repmat(11/12,[1,6]),repmat(3/4,[1,6]),repmat(7/12,[1,6]),repmat(5/12,[1,6]),repmat(1/4,[1,6]),repmat(1/12,[1,6])];
            for i = 1:length(x0)
                s.uMesh          = obj.createBaseDomainIP(x0(i),y0(i));
                obj.perimeter{i} = MinimumThicknessConstraint(s);
            end
        end

        function createVolumeConstraint(obj)
            s.mesh   = obj.mesh;
            s.test = LagrangianFunction.create(obj.mesh,1,'P1');
            s.volumeTarget = 0.5;
            s.uMesh = obj.createBaseDomain();
            v = VolumeConstraint(s);
            obj.volume = v;
        end

        function createCost(obj)
            s.shapeFunctions{1} = obj.compliance;
            s.shapeFunctions{2} = obj.penalty;
            s.weights           = [1,1];
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
            for i = 1:length(obj.perimeter)
                s.shapeFunctions{i+1} = obj.perimeter{i};
            end
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
            s.maxIter        = 3000;
            s.tolerance      = 1e-8;
            s.constraintCase = [{'EQUALITY'},repmat({'INEQUALITY'},[1,36])];
            s.etaNorm        = 0.01;
            s.etaNormMin     = 0.01;
            s.gJFlowRatio    = 4.0;
            s.etaMax         = 50;
            s.etaMaxMin      = 0.1;
            s.primalUpdater  = obj.primalUpdater;
            s.gif            = false;
            s.gifName        = [];
            s.printing       = false;
            s.printName      = [];
            opt = OptimizerNullSpace(s);
            opt.solveProblem();
            obj.optimizer = opt;
        end

        function bc = createBoundaryConditions(obj)
            femReader = FemInputReaderGiD();
            s         = femReader.read(obj.filename);
            sPL       = obj.computeCondition(s.pointload);
            sDir      = obj.computeCondition(s.dirichlet);

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
            s.mesh         = obj.mesh;
            bc = BoundaryConditions(s);
        end
    end

    methods (Static, Access=private)
        function sCond = computeCondition(conditions)
            nodes = @(coor) 1:size(coor,1);
            dirs  = unique(conditions(:,2));
            j     = 0;
            for k = 1:length(dirs)
                rowsDirk = ismember(conditions(:,2),dirs(k));
                u        = unique(conditions(rowsDirk,3));
                for i = 1:length(u)
                    rows   = conditions(:,3)==u(i) & rowsDirk;
                    isCond = @(coor) ismember(nodes(coor),conditions(rows,1));
                    j      = j+1;
                    sCond{j}.domain    = @(coor) isCond(coor);
                    sCond{j}.direction = dirs(k);
                    sCond{j}.value     = u(i);
                end
            end
        end

    end

end