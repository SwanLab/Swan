classdef TopOptPunzonDensity < handle

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
        perimeter
    end

    methods (Access = public)

        function obj = TopOptPunzonDensity()
            obj.init()
            obj.createMesh();
            obj.createDesignVariable();
            obj.createFilter();
            %obj.createPerimeter();
            obj.createMaterialInterpolator();
            obj.createElasticProblem();
            obj.createComplianceFromConstiutive();
            obj.createCompliance();
            obj.createVolumeConstraint();
            obj.createCost();
            obj.createConstraint();
            obj.createPrimalUpdater();
            obj.createOptimizer();

            figure(2)
            saveas(gcf,'Monitoring_PunzonDensity_Intent5.fig');
            obj.designVariable.fun.print('fValues_PunzonDensity_Intent5');
        end

    end

    methods (Access = private)

        function init(obj)
            close all;
        end

        function createMesh(obj)
            file = 'punzon4';
            a.fileName = file;
            s = FemDataContainer(a);
            obj.mesh = s.mesh;
        end

        function createDesignVariable(obj)
            % NON FIXED
            % s.type = 'Full';
            % g      = GeometricalFunction(s);
            % lsFun  = g.computeLevelSetFunction(obj.mesh);
            % s.fun  = lsFun;

            % FIXED
            s.fHandle = @(x) ones(size(x(1,:,:)));
            s.ndimf   = 1;
            s.mesh    = obj.mesh;
            aFun      = AnalyticalFunction(s);
            s.fun     = aFun.project('P1');

            %---------------------------------------------------
            s.mesh = obj.mesh;
            s.type = 'Density';
            s.plotting = true;


            % DESCOMENTAR PER FIXED
            yMin = min(obj.mesh.coord(:,2));
            yMax = max(obj.mesh.coord(:,2));
            zMin     = min(obj.mesh.coord(:,3));
            zMax     = max(obj.mesh.coord(:,3));

            % bolt domain
            c1x= 584.37;    c1y= 28.562;
            c2x= 584.37;    c2y= 68.562;
            hBolt=60;       rBolt=6.5;
            
            bolt1 = @(x) (x(:,1)-c1x).^2 + (x(:,2)-c1y).^2 <= rBolt^2 & ...
                          x(:,3) >= hBolt & x(:,3) <= zMax;
            bolt2 = @(x) (x(:,1)-c2x).^2 + (x(:,2)-c2y).^2 <= rBolt^2 & ...
                          x(:,3) >= hBolt & x(:,3) <= zMax;

            % Guides and bottom fixed
            isBottom = @(x) x(:,3)<= -16; 
            % guide1   = @(x) x(:,2)<= yMin+10.5;
            % guide2   = @(x) x(:,2)>= yMax-10.5;
            guide1   = @(x) x(:,2)<= 20.179;
            guide2   = @(x) x(:,2)>= 76.729;
            s.isFixed  = obj.computeFixedVolumeDomain(...
                         @(x) guide1(x) | guide2(x) | isBottom(x)| bolt1(x) | bolt2(x),...
                         s.type);

            %---------------------------------------------------
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
             s.ub     = 1;
             s.lb     = 0;
             s.tauMax = 1000;
             s.tau    = [];
             obj.primalUpdater = ProjectedGradient(s);
        end

        function createOptimizer(obj)
            s.monitoring     = true;
            s.cost           = obj.cost;
            s.constraint     = obj.constraint;
            s.designVariable = obj.designVariable;
            s.maxIter        = 800;
            s.tolerance      = 1e-8;
            s.constraintCase = {'EQUALITY'};
            s.primalUpdater  = obj.primalUpdater;
            s.etaNorm        = 0.1;
            s.etaNormMin     = 0.02;
            s.gJFlowRatio    = 0.4;
            s.etaMax         = 1;
            s.etaMaxMin      = 0.01;
            s.gif = false;
            s.gifName = [];
            s.printing = true;
            s.printName = 'PunzonDensity';
            opt = OptimizerNullSpace(s);
            opt.solveProblem();
            obj.optimizer = opt;
        end

        function m = createMaterial(obj)
            x = obj.designVariable;
            f = x.obtainDomainFunction();
            f = obj.filter.compute(f{1},1);            
            s.type                 = 'DensityBased';
            %s.fibreOrientation     = '0';
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
end