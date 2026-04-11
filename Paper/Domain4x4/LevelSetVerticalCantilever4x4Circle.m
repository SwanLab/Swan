classdef LevelSetVerticalCantilever4x4Circle < handle

    properties (Access = private)
        mesh
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

        function obj = LevelSetVerticalCantilever4x4Circle(pRelTar)
            obj.init()
            obj.createMesh();
            obj.createDesignVariable();
            obj.createFilter();
            obj.createMaterialInterpolator();
            obj.createElasticProblem();
            obj.createComplianceFromConstiutive();
            obj.createCompliance();
            obj.createPerimeterPenalty();
            obj.createPerimeter(pRelTar);
            obj.createVolumeConstraint();
            obj.createCost();
            obj.createConstraint();
            obj.createPrimalUpdater();
            obj.createOptimizer();

            saveas(gcf,['Paper/Domain4x4/MonitoringLevelSetVerticalCantilever4x4Circle',num2str(pRelTar),'.fig']);
            obj.designVariable.fun.print(['Paper/Domain4x4/LevelSetVerticalCantilever4x4Circle',num2str(pRelTar),'fValues']);
        end

    end

    methods (Access = private)

        function init(obj)
            close all;
        end

        function createMesh(obj)
            obj.mesh = TriangleMesh(1,2,75,150);
        end

        function createDesignVariable(obj)
            s.type             = 'Holes';
            s.dim              = 2;
            s.nHoles           = [50,100];
            s.totalLengths     = [1,2];
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

        function c = createComplianceFromConstiutive(obj)
            s.mesh         = obj.mesh;
            s.stateProblem = obj.physicalProblem;
            c = ComplianceFromConstitutiveTensor(s);
        end

        function createCompliance(obj)
            s.mesh                        = obj.mesh;
            s.filter                      = obj.filter;
            s.complainceFromConstitutive  = obj.createComplianceFromConstiutive();
            s.material                    = obj.createMaterial();
            c = ComplianceFunctional(s);
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
            s.value0      = 6;
            s.signInitial = -0.25;
            s.signFinal   = 0.05;
            s.tarVolume   = 0.4;
            obj.penalty   = InterfaceFunctional(s);
        end

        function uMesh = createBaseDomain(obj)
            levelSet         = -ones(obj.mesh.nnodes,1);
            s.backgroundMesh = obj.mesh;
            s.boundaryMesh   = obj.mesh.createBoundaryMesh();
            uMesh = UnfittedMesh(s);
            uMesh.compute(levelSet);
        end

        function uMesh = createBaseDomainPerimeter(obj,x0,y0)
            s.type             = 'Rectangle';
            s.xCoorCenter      = x0;
            s.yCoorCenter      = y0;
            s.xSide            = 0.25;
            s.ySide            = 0.50;
            g                  = GeometricalFunction(s);
            lsFun              = g.computeLevelSetFunction(obj.mesh);
            sUm.backgroundMesh = obj.mesh;
            sUm.boundaryMesh   = obj.mesh.createBoundaryMesh();
            uMesh              = UnfittedMesh(sUm);
            uMesh.compute(lsFun.fValues);
        end

        function createPerimeter(obj,p)
            sF.mesh       = obj.mesh;
            sF.filterType = 'PDE';
            sF.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
            f             = Filter.create(sF);

            h         = obj.mesh.computeMeanCellSize();
            s.mesh    = obj.mesh;
            s.filter  = f;
            s.epsilon = 3*h;
            s.minEpsilon = 3*h;
            s.value0 = 1;
            s.tarVolume = 0.4;

            tarRef = [0.2372, 0.7385, 0.7385, 0.2373, 0.9665, 0.7976, 0.7974, 0.9665, 0.9162, 1.3799, 1.3797, 0.9159, 0.7075, 0.4713, 0.4714, 0.7077];
            x0 = repmat([0.125,0.375,0.625,0.875],[1,4]);
            y0 = [repmat(1.75,[1,4]),repmat(1.25,[1,4]),repmat(0.75,[1,4]),repmat(0.25,[1,4])];
            for i = 1:length(x0)
                s.uMesh          = obj.createBaseDomainPerimeter(x0(i),y0(i));
                s.target         = p*tarRef(i);
                s.target0        = 100*s.target;
                obj.perimeter{i} = PerimeterConstraint(s);
            end
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
            s.maxIter        = 2500;
            s.tolerance      = 1e-8;
            s.constraintCase = [{'EQUALITY'},repmat({'INEQUALITY'},[1,16])];
            s.etaNorm        = 0.01;
            s.etaNormMin     = 0.01;
            s.gJFlowRatio    = 8.0;
            s.etaMax         = 30;
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
            xMax    = max(obj.mesh.coord(:,1));
            yMax    = max(obj.mesh.coord(:,2));
            isDir   = @(coor)  coor(:,2)==0;
            isForce = @(coor)  coor(:,2)==yMax & coor(:,1)>=0.4*xMax & coor(:,1)<=0.6*xMax;

            sDir{1}.domain    = @(coor) isDir(coor);
            sDir{1}.direction = [1,2];
            sDir{1}.value     = 0;

            sPL{1}.domain    = @(coor) isForce(coor);
            sPL{1}.direction = 1;
            sPL{1}.value     = 1;

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
end