classdef TopOptTestTutorialThermoMechanicalDensity < handle

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
    end

    methods (Access = public)

        function obj = TopOptTestTutorialThermoMechanicalDensity()
            obj.init()
            obj.createMesh();
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
            obj.createOptimizer();
        end

    end

    methods (Access = private)

        function init(obj)
            close all;
        end

        function createMesh(obj)
            x1      = linspace(0,1,60);
            x2      = linspace(0,1,60);
            [xv,yv] = meshgrid(x1,x2);
            [F,V]   = mesh2tri(xv,yv,zeros(size(xv)),'x');
            s.coord  = V(:,1:2);
            s.connec = F;
            obj.mesh = Mesh.create(s);
        end       

        function createDesignVariable(obj)
            s.fHandle = @(x) ones(size(x(1,:,:)));
            s.ndimf   = 1;
            s.mesh    = obj.mesh;
            aFun      = AnalyticalFunction(s);        
            sD.fun      = aFun.project('P1');
            sD.mesh     = obj.mesh;
            sD.type     = 'Density';
            sD.plotting = true;
            dens        = DesignVariable.create(sD);
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
            s.f0   = 90e-3;
            s.f1   = 90;
            s.dim ='2D';
            a = MaterialInterpolator.create(s);
            obj.thermalmaterialInterpolator = a;
         end


        function createMaterialInterpolator(obj)
            E0 =  199.5e3; 
            nu0 = 0.3;
            ndim = obj.mesh.ndim;

            matA.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E0,nu0);
            matA.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E0,nu0,ndim);

            
            E1 =  199.5e9; 
            nu1 = 0.3;
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
            s.alpha = 15.4e-6; 
%             s.source  =  ConstantFunction.create(1,obj.mesh);
%             s.T0 = ConstantFunction.create(25,obj.mesh);   
            s.source  =  ConstantFunction.create(1000,obj.mesh);
            s.T0 = ConstantFunction.create(0,obj.mesh);   
            s.boundaryConditionsThermal = obj.createBoundaryConditionsThermal();
            
           
            fem = ThermoElasticProblem(s);
            obj.physicalProblem = fem;
        end

        function c = createComplianceFromConstiutive(obj)
            s.mesh         = obj.mesh;
            s.stateProblem = obj.physicalProblem;
            s.T0 = ConstantFunction.create(0,obj.mesh);
            c = ComplianceFromConstitutiveTensorThermoElastic(s);
        end

        function createCompliance(obj) 
            s.mesh                        = obj.mesh;
            s.filter                      = obj.filter;
            s.complianceFromConstitutive  = obj.createComplianceFromConstiutive();
            s.material                    = obj.createMaterial();
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
        end

        function createOptimizer(obj)
            s.monitoring     = true;
            s.cost           = obj.cost;
            s.constraint     = obj.constraint;
            s.designVariable = obj.designVariable;
            s.maxIter        = 1000;
            s.tolerance      = 1e-8;
            s.constraintCase = {'INEQUALITY'};
            s.primal         = 'PROJECTED GRADIENT';
            s.etaNorm        = 0.01;
            s.gJFlowRatio    = 1;
            s.gif=false;
            s.gifName='ThermoElastic';
            s.printing=true;
            s.printName='ThermoElastic';
            s.primalUpdater  = obj.primalUpdater;
            opt = OptimizerNullSpace(s);
            opt.solveProblem();
            obj.optimizer = opt;
        end

        function bc = createBoundaryConditionsElastic(obj)
            xMax    = max(obj.mesh.coord(:,1));
            xMin    = min(obj.mesh.coord(:,1));
            yMin    = min(obj.mesh.coord(:,2));
            yMax    = max(obj.mesh.coord(:,2));

%             Cantilever beam
            isDir   = @(coor)  abs(coor(:,1))==0;
            isForce = @(x)  (abs(x(1,:,:))==xMax & abs(x(2,:,:))>=0.4*yMax & abs(x(2,:,:))<=0.6*yMax);

%             Clampled beam
%             isDir   = @(coor)  abs(coor(:,1))== xMin | abs(coor(:,1))==xMax;
%             isForce = @(x)  (abs(x(2,:,:))== yMin & abs(x(1,:,:))>=0.4*xMax & abs(x(1,:,:))<=0.6*xMax);

            sDir{1}.domain    = @(coor) isDir(coor);
            sDir{1}.direction = [1,2];
            sDir{1}.value     = 0;

            [bMesh, ~]  = obj.mesh.createSingleBoundaryMesh();
            sPL{1}.domain = isForce;
            sPL{1}.fun    = ConstantFunction.create([0,-0.5e7/0.2],bMesh);

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

%             isDir   = @(coor) abs(coor(:,2))==yMin & abs(coor(:,1))>=0.4*xMax & abs(coor(:,1))<=0.6*xMax;  
            
            isDir   = @(coor) abs(coor(:,2))==yMin | abs(coor(:,2))==yMax | abs(coor(:,1))==xMin | abs(coor(:,1))==xMax;  

            sDir{1}.domain    = @(coor) isDir(coor);
            sDir{1}.direction = 1;
            sDir{1}.value     = 0;
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
