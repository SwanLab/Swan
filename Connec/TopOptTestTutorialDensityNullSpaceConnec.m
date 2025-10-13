classdef TopOptTestTutorialDensityNullSpaceConnec < handle

    properties (Access = private)
        mesh
        filter
        designVariable
        materialInterpolator
        physicalProblem
        compliance
        volume
        enclosedVoid
        cost
        constraint
        primalUpdater
        optimizer
        perimeter
        filterPDE
    end

    methods (Access = public)

        function obj = TopOptTestTutorialDensityNullSpaceConnec()
            obj.init()
            obj.createMesh();
            obj.createDesignVariable();
            obj.createFilter();
            obj.createMaterialInterpolator();
            obj.createElasticProblem();
            obj.createComplianceFromConstiutive();
            obj.createCompliance();
            obj.createVolumeConstraint();
            obj.createEnclosedVoidFunctionalConstraint();
            obj.createPerimeter();
            obj.createCost();
            obj.createConstraint();
            obj.createPrimalUpdater();
            obj.createOptimizer();
        end

    end

    methods (Access = private)

        function init(obj)
            close all;
        end

        function createMesh(obj)
            %UnitMesh better
            x1      = linspace(0,2,100);
            x2      = linspace(0,1,50);
            [xv,yv] = meshgrid(x1,x2);
            [F,V]   = mesh2tri(xv,yv,zeros(size(xv)),'x');
            s.coord  = V(:,1:2);
            s.connec = F;
            obj.mesh = Mesh.create(s);
        end

        function createDesignVariable(obj)
            % s.fHandle = @(x) ones(size(x(1,:,:)));
            % s.ndimf   = 1;
            % s.mesh    = obj.mesh;
            % aFun      = AnalyticalFunction(s);
            % d = aFun.project('P1');

            s.dim = obj.mesh.ndim;
            s.nHoles = [1 1];
            s.phases = [0 0]; 
            s.phiZero = 0.5;
            s.totalLengths = [1,1];
            s.type         = 'Holes';
            g              = GeometricalFunction(s);
            phiFun         = g.computeLevelSetFunction(obj.mesh);
            ls          = phiFun.fValues;
            lsInclusion = ls;
            sU.backgroundMesh = obj.mesh;
            sU.boundaryMesh   = obj.mesh.createBoundaryMesh;
            uMesh             = UnfittedMesh(sU);
            uMesh.compute(lsInclusion);           
            cFun = CharacteristicFunction.create(uMesh);
            s.filterType = 'LUMP';
            s.mesh  = obj.mesh;
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            f = Filter.create(s);
            d = f.compute(cFun,2);            
            

            % s.fun  = phiFun;
            % s.mesh = obj.mesh;
            % s.type = 'LevelSet';
            % s.plotting = true;
            % ls     = DesignVariable.create(s);
            % obj.designVariable = ls;


            sD.fun      = d;
            sD.mesh     = obj.mesh;
            sD.type     = 'Density';
            sD.plotting = true;
%            sD.pr
            dens        = DesignVariable.create(sD);
            obj.designVariable = dens;
        end

        function createFilter(obj)
            s.filterType = 'LUMP';
            s.mesh  = obj.mesh;
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            f = Filter.create(s);
            obj.filter = f;


            s.filterType        = 'PDE';
            s.boundaryType      = 'Robin';
            s.mesh              = obj.mesh;
            s.trial             = LagrangianFunction.create(obj.mesh,1,'P1');
            f                   = Filter.create(s);
            obj.filterPDE     = f;            
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
            s.filter                      = obj.filter;
            s.complainceFromConstitutive  = obj.createComplianceFromConstiutive();
            s.material                    = obj.createMaterial();
            c = ComplianceFunctional(s);
            obj.compliance = c;
        end

        function uMesh = createBaseDomain(obj)
            levelSet         = -ones(obj.mesh.nnodes,1);
            s.backgroundMesh = obj.mesh;
            s.boundaryMesh   = obj.mesh.createBoundaryMesh();
            uMesh = UnfittedMesh(s);
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

        function createEnclosedVoidFunctionalConstraint(obj)
            %e   = 1e-5;
            h = obj.mesh.computeMeanCellSize;
            e = (h)^2;
            k0  = 1-e;
            k1  = e; 
            m0  = e;
            m1  = 1-e;
            p  = 8;
            s.diffCoef = @(x) k0.*(1-x.^p)+k1*x.^p;
            s.massCoef = @(x) m0.*(1-x.^p)+m1*x.^p;
            s.dm       = @(x) -m0.*p.*x.^(p-1)+m1.*p.*x.^(p-1);
            s.dk       = @(x) -k0.*p.*x.^(p-1)+k1.*p.*x.^(p-1);
            s.mesh   = obj.mesh;
            s.filter = obj.filter;
            s.test   = LagrangianFunction.create(obj.mesh,1,'P1');
            s.uMesh = obj.createBaseDomain();
            v = EnclosedVoidFunctional(s);
            v.computeFunctionAndGradient(obj.designVariable)
            obj.enclosedVoid = v;
        end        

        function createPerimeter(obj)
            eOverhmin     = 10;
            epsilon       = eOverhmin*obj.mesh.computeMeanCellSize();
            s.mesh        = obj.mesh;
            s.filter      = obj.filterPDE;
            s.epsilon     = epsilon;
            s.value0      = 4; % external P
            s.uMesh       = obj.createBaseDomain();
            P             = PerimeterFunctional(s);
            obj.perimeter = P;
        end   
  

        function createCost(obj)
            s.shapeFunctions{1} = obj.enclosedVoid;
            s.shapeFunctions{2} = obj.perimeter;
            s.weights           = [1 100];
            s.Msmooth           = obj.createMassMatrix();
            obj.cost            = Cost(s);
        end

        function M = createMassMatrix(obj)
            test  = LagrangianFunction.create(obj.mesh,1,'P1');
            trial = LagrangianFunction.create(obj.mesh,1,'P1'); 
            M = IntegrateLHS(@(u,v) DP(v,u),test,trial,obj.mesh,'Domain');
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

            % 
            % s.mesh = obj.mesh;
            % obj.primalUpdater = SLERP(s);            
        end

        function createOptimizer(obj)
            % s.monitoring     = true;
            % s.cost           = obj.cost;
            % s.constraint     = obj.constraint;
            % s.designVariable = obj.designVariable;
            % s.maxIter        = 3000;
            % s.tolerance      = 1e-8;
            % s.constraintCase = {'EQUALITY'};
            % s.primalUpdater  = obj.primalUpdater;
            % s.etaNorm        = 0.02;
            % s.etaNormMin     = 0.02;
            % s.gJFlowRatio    = 0.1;
            % s.etaMax         = 1;
            % s.etaMaxMin      = 0.01;
            % opt = OptimizerNullSpace(s);
            % opt.solveProblem();
            % obj.optimizer = opt;


            s.nConstraints   = 1;
            l                = DualVariable(s);            


            s.monitoring     = true;
            s.cost           = obj.cost;
            s.constraint     = obj.constraint;
            s.designVariable = obj.designVariable;
            s.dualVariable   = l;
            s.maxIter        = 500;
            s.tolerance      = 1e-8;
            s.constraintCase = {'EQUALITY'};
            s.ub             = 1;
            s.lb             = 0;
            s.volumeTarget   = 0.4;
            s.primal         = 'PROJECTED GRADIENT';
            opt              = OptimizerMMA(s);
            opt.solveProblem();
            obj.optimizer = opt;

            % 
            % s.monitoring     = true;
            % s.cost           = obj.cost;
            % s.constraint     = obj.constraint;
            % s.designVariable = obj.designVariable;
            % s.dualVariable   = l;
            % s.maxIter        = 1000;
            % s.tolerance      = 1e-8;
            % s.constraintCase = {'EQUALITY'};
            % s.tauMax = 1000;
            % s.primal         = 'PROJECTED GRADIENT';
            % s.ub             = 1;
            % s.lb             = 0;
            % s.rho            = obj.designVariable;
            % s.volumeTarget   = 0.4;
            % opt = OptimizerAugmentedLagrangian(s);
            % opt.solveProblem();
            % obj.optimizer = opt;
            % opt.solveProblem();
            % obj.optimizer = opt;
        end

        function bc = createBoundaryConditions(obj)
            xMax    = max(obj.mesh.coord(:,1));
            yMax    = max(obj.mesh.coord(:,2));
            isDir   = @(coor)  abs(coor(:,1))==0;
            isForce = @(coor)  (abs(coor(:,1))==xMax & abs(coor(:,2))>=0.4*yMax & abs(coor(:,2))<=0.6*yMax);

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
            % for i = 1:numel(sPL)
            %     pl = PointLoad(obj.mesh, sPL{i});
            %     pointloadFun = [pointloadFun, pl];
            % end
            s.pointloadFun = pointloadFun;

            s.periodicFun  = [];
            s.mesh         = obj.mesh;
            bc = BoundaryConditions(s);
        end
    end
end
