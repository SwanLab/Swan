classdef TopOptTestTutorialThermoMechanicalBatteryLevelSet < handle

    properties (Access = private)
        mesh
        filter
        designVariable
        materialInterpolator
        thermalmaterialInterpolator 
        physicalProblem
        compliance
        volume
        cost
        constraint
        primalUpdater
        optimizer
        chiB    %Battery
        chiAl  %Corona
        chiB0
        kappaB
        dualVariable
    end

    methods (Access = public)

        function obj = TopOptTestTutorialThermoMechanicalBatteryLevelSet()
            obj.init()
            obj.createMesh();
            obj.createBatteryDomain();
            obj.createNonDesignableDomain();
            obj.createDesignVariable();
            obj.createFilter();
            obj.createMaterialInterpolator();
            obj.createThermalMaterialInterpolator(); 
            obj.createThermoElasticProblem();
            obj.createComplianceFromConstiutive();
            obj.createCompliance();
            obj.createVolumeConstraint();
            obj.createCost();
            obj.createConstraint();
             obj.createPrimalUpdater();
            %obj.createDualVariable();
            obj.createOptimizer();
        end

    end

    methods (Access = private)

        function init(obj)
            close all;
        end

        function createMesh(obj)
            x1      = linspace(0,1,100);
            x2      = linspace(0,1,100);
            [xv,yv] = meshgrid(x1,x2);
            [F,V]   = mesh2tri(xv,yv,zeros(size(xv)),'x');
            s.coord  = V(:,1:2);
            s.connec = F;
            obj.mesh = Mesh.create(s);
        end

         function createBatteryDomain(obj)
             s.type        = 'Circles';
             s.r           = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01].*10;
             s.x0 = [0, 0, 0, 0.04, 0.08, 0.04, 0.08, 0.04, 0.08].*10;
             s.y0 = [0, 0.04, 0.08, 0, 0, 0.04, 0.08, 0.08, 0.04].*10;
             g                = GeometricalFunction(s);
             phiFun           = g.computeLevelSetFunction(obj.mesh);
             phi              = phiFun.fValues;

             sm.backgroundMesh = obj.mesh;
             sm.boundaryMesh   = obj.mesh.createBoundaryMesh;
             uMesh              = UnfittedMesh(sm);
             uMesh.compute(phi);

             obj.chiB     = CharacteristicFunction.create(uMesh).project('P1');
         end        

        function createNonDesignableDomain(obj)
             %corona di alluminio
             sR.type        = 'Circles';
             sR.r           = [0.013, 0.013, 0.013, 0.013, 0.013, 0.013, 0.013, 0.013, 0.013].*10;
             sR.x0 = [0, 0, 0, 0.04, 0.08, 0.04, 0.08, 0.04, 0.08].*10;
             sR.y0 = [0, 0.04, 0.08, 0, 0, 0.04, 0.08, 0.08, 0.04].*10;
             gR                = GeometricalFunction(sR);
             phiFunR           = gR.computeLevelSetFunction(obj.mesh);
             phiR              = phiFunR.fValues;

             smR.backgroundMesh = obj.mesh;
             smR.boundaryMesh   = obj.mesh.createBoundaryMesh;
             uMeshR              = UnfittedMesh(smR);
             uMeshR.compute(phiR);

             obj.chiAl     = CharacteristicFunction.create(uMeshR).project('P1'); %- obj.chiB;


        end

        function createDesignVariable(obj)
            s.type = 'Full';
            g      = GeometricalFunction(s);
            lsFun  = g.computeLevelSetFunction(obj.mesh);
            s.fun  = lsFun;
            s.mesh = obj.mesh;
            s.type = 'LevelSet';
            s.plotting = true;
%             obj.chiB0 = project(obj.chiB,'P1');
            plot(obj.chiB);
            
            x0 = [0, 0, 0, 0.04, 0.08, 0.04, 0.08, 0.04, 0.08].*10;
            y0 = [0, 0.04, 0.08, 0, 0, 0.04, 0.08, 0.08, 0.04].*10;
            centers = [x0;y0]'; r = 0.13;
            isNonDesign =  @(coor) any( (coor(:,1) - centers(:,1)').^2 + ...
                                       (coor(:,2) - centers(:,2)').^2 <= r^2, 2)| coor(:,1) > 0.97 | coor(:,2) > 0.97;

            s.isFixed.nodes  =  isNonDesign(obj.mesh.coord);
            s.isFixed.values = -1; 

            dens        = DesignVariable.create(s);
            obj.designVariable = dens;
        end

        function createFilter(obj)
            s.filterType = 'LUMP';
            s.mesh  = obj.mesh;
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            f = Filter.create(s);
            obj.filter = f;
        end

         function createThermalMaterialInterpolator(obj) % Conductivity
            s.interpolation  = 'SimpAllThermal';   
            s.f0   = 1e-2; %ConstantFunction.create(1e-3,obj.mesh); %1e-2;
            s.f1   = 1.0; %ConstantFunction.create(1,obj.mesh); 
            s.dim ='2D';
            a = MaterialInterpolator.create(s);
            obj.thermalmaterialInterpolator = a;
         end

        function createMaterialInterpolator(obj)
            
            E0 =  ConstantFunction.create(1e-3,obj.mesh);
            nu0 = 1/3;
            ndim = obj.mesh.ndim;

            matA.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E0,nu0);
            matA.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E0,nu0,ndim);

            Ea = ConstantFunction.create(1,obj.mesh); % Aluminium
            Eb = ConstantFunction.create(1.5/68,obj.mesh); % Battery
            
            obj.createBatteryDomain();
            E1 = Ea.*(1 - obj.chiB) + Eb.*obj.chiB; 

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
            x = obj.designVariable;
            f = x.obtainDomainFunction();
            f = obj.filter.compute(f{1},1);
            s.type                 = 'DensityBased';
            s.density              = f;
            s.materialInterpolator = obj.materialInterpolator;
            s.dim                  = '2D';
            s.mesh                 = obj.mesh;
            m = Material.create(s);
        end

       function createThermoElasticProblem(obj)
            s.mesh = obj.mesh;
            s.dim = '2D';
            s.interpolationType = 'LINEAR';
            s.solverType = 'REDUCED';
            s.solverMode = 'DISP';
            s.solverCase = DirectSolver();

            % Elastic
            s.scale = 'MACRO';
            s.material = obj.createMaterial();
            s.boundaryConditionsElastic = obj.createBoundaryConditionsElastic();

            % Thermal
            s.materialInterpolator = obj.thermalmaterialInterpolator;
            s.alpha = 4.0; 
            s.source  =  ConstantFunction.create(1,obj.mesh).*obj.chiB.project('P1'); %P=2.5
            s.T0 = ConstantFunction.create(0,obj.mesh);   
            s.boundaryConditionsThermal = obj.createBoundaryConditionsThermal();
            
           
            fem = ThermoElasticProblem(s);
            obj.physicalProblem = fem;
        end

        function c = createComplianceFromConstiutive(obj)
            s.T0 = ConstantFunction.create(0,obj.mesh);
            s.mesh         = obj.mesh;
            s.stateProblem = obj.physicalProblem;
            c = ComplianceFromConstitutiveTensorThermoElastic(s);
        end

        function createCompliance(obj) 
            s.mesh                        = obj.mesh;
            s.filter                      = obj.filter;
            s.complianceFromConstitutive  = obj.createComplianceFromConstiutive();
            s.material                    = obj.createMaterial();
            s.chiB                        = obj.chiB;
            s.kappaB                      = 1.25/220;
            s.conductivity                = obj.thermalmaterialInterpolator();
            c = ComplianceFunctionalThermoElastic(s);
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
            s.volumeTarget = 0.5;
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
            s.maxIter        = 400;
            s.tolerance      = 1e-8;
            s.constraintCase = {'EQUALITY'};
            s.primalUpdater  = obj.primalUpdater;
            s.etaNorm        = 0.02;
            s.etaNormMin     = 0.02;
            s.gJFlowRatio    = 1;
            s.etaMax         = 1;
            s.etaMaxMin      = 0.01;
            s.gif            =true;
            s.gifName        ='BatteryLevelSet';
            s.printing       =true;
            s.printName      ='Battery Level Set';
            opt = OptimizerNullSpace(s);
            opt.solveProblem();
            obj.optimizer = opt;
        end

        function bc = createBoundaryConditionsElastic(obj)
            xMin    = min(obj.mesh.coord(:,1));
            xMax    = max(obj.mesh.coord(:,1));
            yMin    = min(obj.mesh.coord(:,2));
            yMax    = max(obj.mesh.coord(:,2));

            % using a 1/4 of the structure exploiting the simmetries
            isDirDown    = @(coor)  abs(coor(:,2))==yMin;     
            isDirLeft    = @(coor)  abs(coor(:,1))==xMin;   
            isForceRight = @(x)  abs(x(1,:,:))==xMax;  % up and right Neumann!
            isForceUp    = @(x)  abs(x(2,:,:))==yMax;
           
            sDir{1}.domain    = @(coor) isDirLeft(coor);
            sDir{1}.direction = 1;
            sDir{1}.value     = 0;
            sDir{2}.domain    = @(coor) isDirDown(coor);
            sDir{2}.direction = 2;
            sDir{2}.value     = 0;

            [bMesh, ~]  = obj.mesh.createSingleBoundaryMesh();
            sPL{1}.domain = isForceRight;
            sPL{1}.fun    = ConstantFunction.create([-1,0],bMesh);
            sPL{2}.domain = isForceUp;
            sPL{2}.fun    = ConstantFunction.create([0,-1],bMesh);

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

            s.periodicFun  = [];
            s.mesh         = obj.mesh;
            bc = BoundaryConditions(s);
        end

        function bcT = createBoundaryConditionsThermal(obj) 
            xMin    = min(obj.mesh.coord(:,1));
            yMin    = min(obj.mesh.coord(:,2));
            xMax    = max(obj.mesh.coord(:,1));
            yMax    = max(obj.mesh.coord(:,2));
            
            % dirichlet on up and right
            isDir   = @(coor)  abs(coor(:,2))==yMax | abs(coor(:,1))==xMax;  

            sDir{1}.domain    = @(coor) isDir(coor);
            sDir{1}.direction = 1;
            sDir{1}.value     = 0; %20;
            sDir{1}.ndim = 1;
            
            dirichletFun = [];
            for i = 1:numel(sDir)
                dir = DirichletCondition(obj.mesh, sDir{i});
                dirichletFun = [dirichletFun, dir];
            end
            s.dirichletFun = dirichletFun;
            s.pointloadFun = [];

            s.periodicFun  = [];
            s.mesh         = obj.mesh;
            bcT = BoundaryConditions(s); 
         end
    end
end
