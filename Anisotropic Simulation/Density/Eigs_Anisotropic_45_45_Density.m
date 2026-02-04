classdef Eigs_Anisotropic_45_45_Density < handle

    properties (Access = private)
        mesh
        filter
        filterConnect
        designVariable
        materialInterpolator
        massInterpolator
        physicalProblem
        compliance
        volume
        cost
        constraint
        primalUpdater
        optimizer
        perimeter
        minimumEigenValue
        lambda1min

    end

    methods (Access = public)

        function obj = Eigs_Anisotropic_45_45_Density()
            obj.lambda1min = 1e-5; % To change with the 50 Hz minimum
            obj.init()
            obj.createMesh();
            obj.createDesignVariable();
            obj.createFilter();
            obj.createFilterConnectivity();
            %obj.createPerimeter();
            obj.createMaterialInterpolator();
            obj.createMassInterpolator();
            obj.createElasticProblem();
            obj.createComplianceFromConstiutive();
            obj.createCompliance();
            obj.createVolumeConstraint();
            obj.createEigenValueConstraint();  
            obj.createCost();
            obj.createConstraint();
            obj.createPrimalUpdater();
            obj.createOptimizer();

            % Save monitoring and desginVariable fValues
            figure(2);
            saveas(gcf,'Monitoring_45_45_Density_Eigs_Lambda_1minus5.fig');
            obj.designVariable.fun.print('fValues_45_45_Density_Eigs_Lambda_1minus5');
        end

    end

    methods (Access = private)

        function init(obj)
            close all;
        end

        function createMesh(obj)
            %UnitMesh better
            % Cantilever beam
            x1      = linspace(0,2,150);
            x2      = linspace(0,1,75);
            % MBB Beam
            % x1      = linspace(0,6,500);
            % x2      = linspace(0,1,75);
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

        function createFilterConnectivity(obj)
            s.filterType = 'LUMP';
            s.mesh       = obj.mesh;
            s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
            f            = Filter.create(s);
            obj.filterConnect = f; 
        end
  
        function createMaterialInterpolator(obj)
            type = '-45_45';
            s.C1 = Cvoigt.create(type);
            s.C0 = s.C1*1e-3; % This is not necessary
            s.mesh = obj.mesh;
            s.interpolation  = 'SIMP_P3_ANISOTROPIC';
            s.dim            = '2D';

            m = MaterialInterpolator.create(s);
            obj.materialInterpolator = m;
        end


        function createMassInterpolator(obj)
            s.interpolation  = 'SIMPThermal';                              
            s.f0   = 1e-3;
            s.f1   = 1;
            s.pExp = 1;
            a = MaterialInterpolator.create(s);
            obj.massInterpolator = a;            
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

         function createEigenValueConstraint(obj)                           
            s.mesh              = obj.mesh;
            s.designVariable    = obj.designVariable;
            s.filter            = obj.filterConnect;
            s.boundaryConditions = obj.createEigenvalueBoundaryConditions();
            s.material           = obj.createMaterial(); 
            s.massInterpolator   = obj.massInterpolator; 
            s.targetEigenValue  = obj.lambda1min;    
            s.isCompl           = true;
            obj.minimumEigenValue = StiffnessEigenModesConstraint(s);
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
            s.shapeFunctions{2} = obj.minimumEigenValue; 
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
            s.maxIter        = 600;
            s.tolerance      = 1e-8;
            s.constraintCase = {'EQUALITY', 'INEQUALITY'};
            s.primalUpdater  = obj.primalUpdater;
            s.etaNorm        = 0.02;
            s.etaNormMin     = 0.02;
            s.gJFlowRatio    = 1;
            s.etaMax         = 1;
            s.etaMaxMin      = 0.01;
            %s.type           = '0';
            s.gif = true;
            s.gifName = 'Eigs_Gif_45_45_Density_Lambda_1minus5';
            s.printing = false;
            s.printName = 'Eigs_Results_45_45_Density_Lambda_1minus5';
            opt = OptimizerNullSpace(s);
            opt.solveProblem();
            obj.optimizer = opt;
        end

        function m = createMaterial(obj)
            x = obj.designVariable;
            f = x.obtainDomainFunction();
            f = obj.filter.compute(f{1},1);            
            s.type                 = 'DensityBasedMaterialAnisotropic';
            %s.fibreOrientation     = '0';
            s.density              = f;
            s.materialInterpolator = obj.materialInterpolator;
            s.dim                  = '2D';
            s.mesh                 = obj.mesh;
            m = Material.create(s);
        end

        function bc = createBoundaryConditions(obj)
            % Cantilever beam
            xMax    = max(obj.mesh.coord(:,1));
            yMax    = max(obj.mesh.coord(:,2));
            isDir   = @(coor)  abs(coor(:,1))==0;
            isForce = @(coor)  (abs(coor(:,1))==xMax & abs(coor(:,2))>=0.4*yMax & abs(coor(:,2))<=0.6*yMax);
            
            % MBB beam
            % xMax    = max(obj.mesh.coord(:,1));
            % yMax    = max(obj.mesh.coord(:,2));
            % isDir   = @(coor)  (abs(coor(:,1))==0 & abs(coor(:,2)) == 0);
            % isDir2  = @(coor)  (abs(coor(:,1))==xMax & abs(coor(:,2)) == 0);
            % isForce = @(coor)  (abs(coor(:,2))==yMax & abs(coor(:,1))>=0.4*xMax & abs(coor(:,1))<=0.6*xMax);


            sDir{1}.domain    = @(coor) isDir(coor);
            sDir{1}.direction = [1,2]; % Cantilever--> [1,2]   MBB--> 2
            sDir{1}.value     = 0;
            
            % Comentar sDir 2 quan es faci cantilever beam
            % sDir{2}.domain    = @(coor) isDir2(coor);
            % sDir{2}.direction = [1,2];
            % sDir{2}.value     = 0;

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
                pl = TractionLoad(obj.mesh, sPL{i},'DIRAC');
                pointloadFun = [pointloadFun, pl];
            end
            s.pointloadFun = pointloadFun;

            s.periodicFun  = [];
            s.mesh = obj.mesh;
            bc = BoundaryConditions(s);
        end

         function  bc = createEigenvalueBoundaryConditions(obj)
            xMin    = min(obj.mesh.coord(:,1));
            yMin    = min(obj.mesh.coord(:,2));
            xMax    = max(obj.mesh.coord(:,1));
            yMax    = max(obj.mesh.coord(:,2));
            isLeft  = @(coor) abs(coor(:,1))==xMin;
            isRight = @(coor) abs(coor(:,1))==xMax;
            isFront = @(coor) abs(coor(:,2))==yMin;
            isBack = @(coor) abs(coor(:,2))== yMax;
            isDir   = @(coor) isLeft(coor);% | isRight(coor) | isFront(coor) | isBack(coor);  
            sDir{1}.domain    = @(coor) isDir(coor);
            sDir{1}.direction = [1,2];
            sDir{1}.value     = 0;
            sDir{1}.ndim = 2;
            
            dirichletFun = [];
            for i = 1:numel(sDir)
                dir = DirichletCondition(obj.mesh, sDir{i});
                dirichletFun = [dirichletFun, dir];
            end
            s.dirichletFun = dirichletFun;
            s.pointloadFun = [];

            s.periodicFun  = [];
            s.mesh         = obj.mesh;
            bc = BoundaryConditions(s);  
         end

         
    end
end