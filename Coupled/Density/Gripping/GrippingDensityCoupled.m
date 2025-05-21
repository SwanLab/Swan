classdef GrippingDensityCoupled < handle

    properties (Access = private)
        etaSt

        mesh
        filename
        filterCompliance
        designVariable
        materialInterpolator
        physicalProblem
        compliance
        globalPer
        volume
        perimeter
        cost
        constraint
        primalUpdater
        optimizer
    end

    methods (Access = public)

        function obj = GrippingDensityCoupled(eta)
            obj.init(eta)
            obj.createMesh();
            obj.createDesignVariable();
            obj.createFilterCompliance();
            obj.createMaterialInterpolator();
            obj.createElasticProblem();
            obj.createNonSelfAdjCompliance();
            obj.createGlobalPerimeter();
            obj.createVolumeConstraint();
            obj.createPerimeter();
            obj.createCost();
            obj.createConstraint();
            obj.createPrimalUpdater();
            obj.createOptimizer();

            d = obj.designVariable;
            info = ['Eta',num2str(obj.etaSt)];
            saveas(gcf,['Coupled/Density/Gripping/Monit',info,'.fig']);
            save(['Coupled/Density/Gripping/DesVar',info,'.mat'],'d');
        end

    end

    methods (Access = private)

        function init(obj,eta)
            close all;
            obj.etaSt = eta;
        end

        function createMesh(obj)
            file = 'Gripping';
            obj.filename = file;
            a.fileName = file;
            s = FemDataContainer(a);
            obj.mesh = s.mesh;        
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

        function f = createFilterPerimeter(obj)
            s.trial   = LagrangianFunction.create(obj.mesh,1,'P1');
            s.mesh    = obj.mesh;
            s.LHStype = 'StiffnessMass';
            f         = FilterPDE(s);
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

        function createNonSelfAdjCompliance(obj)
            s.mesh         = obj.mesh;
            s.filter       = obj.filterCompliance;
            s.material     = obj.createMaterial();
            s.stateProblem = obj.physicalProblem;
            s.filename     = obj.filename;
            c = NonSelfAdjointComplianceFunctional(s);
            obj.compliance = c;
        end

        function createGlobalPerimeter(obj)
            PMax     = 9*0.9;
            s.target = PMax;

            s.mesh        = obj.mesh;
            s.epsilon     = 6*obj.mesh.computeMeanCellSize();
            s.minEpsilon  = 1.5*obj.mesh.computeMeanCellSize();
            s.value0      = 4;
            s.uMesh       = obj.createBaseGlobalDomain();
            s.filter      = obj.createFilterPerimeter();
            obj.globalPer = PerimeterConstraint(s);
        end

        function createVolumeConstraint(obj)
            s.mesh   = obj.mesh;
            s.filter = obj.filterCompliance;
            s.gradientTest = LagrangianFunction.create(obj.mesh,1,'P1');
            s.volumeTarget = 0.35;
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

        function base = createBaseDomain(obj,x0,y0)
            s.type             = 'Circle';
            s.xCoorCenter      = x0;
            s.yCoorCenter      = y0;
            s.radius           = 0.02;
            g                  = GeometricalFunction(s);
            lsFun              = g.computeLevelSetFunction(obj.mesh);
            sUm.backgroundMesh = obj.mesh;
            sUm.boundaryMesh   = obj.mesh.createBoundaryMesh();
            base               = UnfittedMesh(sUm);
            base.compute(lsFun.fValues);
        end

        function createPerimeter(obj)
            s.mesh       = obj.mesh;
            s.epsilon    = 6*obj.mesh.computeMeanCellSize();
            s.minEpsilon = 1.5*obj.mesh.computeMeanCellSize();
            s.value0     = 2*pi*0.02;

            PMax     = 9*0.9;
            L        = sqrt(2);
            l        = 0.04;
            l0       = 0.05; % Minimum length scale
            s.target = max(0,(PMax/L^2)*(l^2-l0^2));

            s.uMesh  = obj.createBaseDomain(0.28,0.5);
            s.filter = obj.createFilterPerimeter();
            P{1}     = PerimeterConstraint(s);

            s.uMesh  = obj.createBaseDomain(0.475,0.745);
            s.filter = obj.createFilterPerimeter();
            P{2}     = PerimeterConstraint(s);

            s.uMesh  = obj.createBaseDomain(0.475,0.255);
            s.filter = obj.createFilterPerimeter();
            P{3}     = PerimeterConstraint(s);

            s.uMesh  = obj.createBaseDomain(0.75,0.5);
            s.filter = obj.createFilterPerimeter();
            P{4}     = PerimeterConstraint(s);

            obj.perimeter = P;
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
            s.shapeFunctions{2} = obj.globalPer;
            s.shapeFunctions{3} = obj.perimeter{1};
            s.shapeFunctions{4} = obj.perimeter{2};
            s.shapeFunctions{5} = obj.perimeter{3};
            s.shapeFunctions{6} = obj.perimeter{4};
            s.Msmooth      = obj.createMassMatrix();
            obj.constraint = Constraint(s);
        end

        function createPrimalUpdater(obj)
            s.ub     = 1;
            s.lb     = 0;
            s.tauMax = 1000;
            obj.primalUpdater = ProjectedGradient(s);
        end

        function createOptimizer(obj)
            s.monitoring     = true;
            s.cost           = obj.cost;
            s.constraint     = obj.constraint;
            s.designVariable = obj.designVariable;
            s.maxIter        = 3000;
            s.tolerance      = 1e-8;
            s.constraintCase = repmat({'INEQUALITY'},[1,6]);
            s.primal         = 'PROJECTED GRADIENT';
            s.etaNorm        = 0.02;
            s.gJFlowRatio    = obj.etaSt;
            s.primalUpdater  = obj.primalUpdater;
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
                pl = PointLoad(obj.mesh, sPL{i});
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
