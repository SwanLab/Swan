classdef CantileverDensityCoupled < handle

    properties (Access = private)
        etaSt

        mesh
        filterCompliance
        filterGlobPer
        filterSegPer
        designVariable
        materialInterpolator
        physicalProblem
        compliance
        globalPer
        volume
        segPerimeter
        cost
        constraint
        primalUpdater
        optimizer
    end

    methods (Access = public)

        function obj = CantileverDensityCoupled(eta)
            obj.init(eta)
            obj.createMesh();
            obj.createDesignVariable();
            obj.createFilterCompliance();
            obj.createFilterGlobalPerimeter();
            obj.createFilterSegmentPerimeter();
            obj.createMaterialInterpolator();
            obj.createElasticProblem();
            obj.createCompliance();
            obj.createGlobalPerimeter();
            obj.createVolumeConstraint();
            obj.createSegmentPerimeter();
            obj.createCost();
            obj.createConstraint();
            obj.createPrimalUpdater();
            obj.createOptimizer();

            d = obj.designVariable;
            info = ['Eta',num2str(obj.etaSt)];
            saveas(gcf,['Coupled/Density/Cantilever/Monit',info,'.fig']);
            save(['Coupled/Density/Cantilever/DesVar',info,'.mat'],'d');
        end

    end

    methods (Access = private)

        function init(obj,eta)
            close all;
            obj.etaSt = eta;
        end

        function createMesh(obj)
            obj.mesh = TriangleMesh(2,1,220,110);
        end

        function createDesignVariable(obj)
            s.fHandle = @(x) ones(size(x(1,:,:)));
            s.ndimf   = 1;
            s.mesh    = obj.mesh;
            aFun      = AnalyticalFunction(s);
            
            sD.fun      = aFun.project('P1');
            sD.mesh     = obj.mesh;
            sD.type     = 'Density';
            sD.plotting = false;
            dens        = DesignVariable.create(sD);
            obj.designVariable = dens;
        end

        function createFilterCompliance(obj)
            s.filterType = 'LUMP';
            s.mesh  = obj.mesh;
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            f = Filter.create(s);
            obj.filterCompliance = f;
        end

        function createFilterGlobalPerimeter(obj)
            s.trial   = LagrangianFunction.create(obj.mesh,1,'P1');
            s.mesh    = obj.mesh;
            s.LHStype = 'StiffnessMass';
            obj.filterGlobPer = FilterPDE(s);
        end

        function createFilterSegmentPerimeter(obj)
            s.mesh  = obj.mesh;
            s.theta = 90;
            s.alpha = 8;
            s.beta  = 0;
            obj.filterSegPer = NonLinearFilterSegment(s);
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
            s.solverCase = 'DIRECT';
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
            s.filter                      = obj.filterCompliance;
            s.complainceFromConstitutive  = obj.createComplianceFromConstiutive();
            s.material                    = obj.createMaterial();
            c = ComplianceFunctional(s);
            obj.compliance = c;
        end

        function createGlobalPerimeter(obj)
            s.mesh        = obj.mesh;
            s.epsilon     = 4*obj.mesh.computeMeanCellSize();
            s.value0      = 6;
            s.uMesh       = obj.createBaseGlobalDomain();
            s.filter      = obj.filterGlobPer;
            obj.globalPer = PerimeterFunctional(s);
        end

        function createVolumeConstraint(obj)
            s.mesh   = obj.mesh;
            s.filter = obj.filterCompliance;
            s.gradientTest = LagrangianFunction.create(obj.mesh,1,'P1');
            s.volumeTarget = 0.4;
            v = VolumeConstraint(s);
            obj.volume = v;
        end

        function base = createBaseGlobalDomain(obj)
            s.type             = 'Full';
            g                  = GeometricalFunction(s);
            lsFun              = g.computeLevelSetFunction(obj.mesh);
            sUm.backgroundMesh = obj.mesh;
            sUm.boundaryMesh   = obj.mesh.createBoundaryMesh();
            base               = UnfittedMesh(sUm);
            base.compute(lsFun.fValues);
        end

        function createSegmentPerimeter(obj)
            s.mesh       = obj.mesh;
            s.epsilon    = obj.mesh.computeMeanCellSize();
            s.minEpsilon = obj.mesh.computeMeanCellSize();
            s.value0     = 2;
            s.target     = 1000; % esto sera un cero numerico
            s.uMesh      = obj.createBaseGlobalDomain();
            s.filter     = obj.filterSegPer;
            obj.segPerimeter = PerimeterConstraint(s);
        end

        function createCost(obj)
            s.shapeFunctions{1} = obj.compliance;
            s.shapeFunctions{2} = obj.globalPer;
            s.weights           = [1,0.80]; % 0.25 hay más barras que tmb cumplen min length, 0.50 + 2 barritas pequeñas
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
            s.shapeFunctions{2} = obj.segPerimeter;
            s.Msmooth      = obj.createMassMatrix();
            obj.constraint = Constraint(s);
        end

        function createPrimalUpdater(obj)
            s.ub     = 1;
            s.lb     = 0;
            s.tauMax = 1000;
            s.tau = [];
            obj.primalUpdater = ProjectedGradient(s);
        end

        function createOptimizer(obj)
            s.monitoring     = true;
            s.cost           = obj.cost;
            s.constraint     = obj.constraint;
            s.designVariable = obj.designVariable;
            s.maxIter        = 3000;
            s.tolerance      = 1e-8;
            s.constraintCase = repmat({'INEQUALITY'},[1,2]);
            s.primal         = 'PROJECTED GRADIENT';
            s.etaNorm        = 0.02;
            s.gJFlowRatio    = obj.etaSt;
            s.primalUpdater  = obj.primalUpdater;
            opt = OptimizerNullSpace(s);
            opt.solveProblem();
            obj.optimizer = opt;
        end

        function bc = createBoundaryConditions(obj)
            xMax    = max(obj.mesh.coord(:,1));
            yMax    = max(obj.mesh.coord(:,2));
            isDir   = @(coor)  coor(:,1)==0;
            isForce = @(coor)  coor(:,1)==xMax & coor(:,2)>=0.4*yMax & coor(:,2)<=0.6*yMax;

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
                pl = PointLoad(obj.mesh, sPL{i});
                pointloadFun = [pointloadFun, pl];
            end
            s.pointloadFun = pointloadFun;

            s.periodicFun  = [];
            s.mesh         = obj.mesh;
            bc = BoundaryConditions(s);
        end
    end

end
