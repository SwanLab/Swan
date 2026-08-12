classdef TutorialFirstMMA < handle

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
        dualVariable
        optimizer
    end

    methods (Access = public)

        function obj = TutorialFirstMMA()
            obj.init()
            obj.createMesh();
            obj.createDesignVariable();
            obj.createFilter();
            % obj.createFilterPerimeter();
            %obj.createMaterialInterpolator();
            obj.createElasticProblem();
            obj.createBaseDomain();
            obj.createComplianceFromConstiutive();
            obj.createCompliance();
            % obj.createVolumeFunctional();
            % obj.createComplianceConstraint();
            obj.createVolumeConstraint();
            % obj.createPerimeter();
            obj.createCost();
            obj.createConstraint();
            obj.createDualVariable();
            obj.createOptimizer();
        end

    end

    methods (Access = private)

        function init(obj)
            close all;
        end

        function createMesh(obj)
           obj.mesh = TriangleMesh(2,1,100,50);
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
            rho        = DesignVariable.create(sD);

            obj.designVariable = rho;
        end

        function createFilter(obj)
            s.filterType = 'LUMP';
            s.mesh  = obj.mesh;
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            f = Filter.create(s);
            obj.filter = f;
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

        function uMesh = createBaseDomain(obj)
            levelSet         = -ones(obj.mesh.nnodes,1);
            s.backgroundMesh = obj.mesh;
            s.boundaryMesh   = obj.mesh.createBoundaryMesh();
            uMesh = UnfittedMesh(s);
            uMesh.compute(levelSet);
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
        



        % function createComplianceConstraint(obj)
        %     s.mesh                       = obj.mesh;
        %     s.filter                     = obj.filterCompliance;
        %     s.complainceFromConstitutive = obj.createComplianceFromConstiutive();
        %     s.material                   = obj.createMaterial();
        %     s.complianceTarget           = 3;
        %     c = ComplianceConstraint(s);
        %     obj.compliance = c;
        % end
        % 
        function createVolumeConstraint(obj)
            s.mesh   = obj.mesh;
            s.filter = obj.filter;
            s.test = LagrangianFunction.create(obj.mesh,1,'P1');
            s.volumeTarget = 0.4;
            s.uMesh = obj.createBaseDomain();
            v = VolumeConstraint(s);
            obj.volume = v;
        end

        % function createPerimeter(obj)
        %     eOverhmin     = 10; % 10
        %     epsilon       = eOverhmin*obj.mesh.computeMeanCellSize();
        %     s.mesh        = obj.mesh;
        %     s.filter      = obj.filterPerimeter;
        %     s.epsilon     = epsilon;
        %     s.value0      = 6; % external Perimeter
        %     P             = PerimeterFunctional(s);
        %     obj.perimeter = P;
        % end

        function createCost(obj)
            s.shapeFunctions{1} = obj.compliance;
            s.weights           = 1;
            s.Msmooth           = obj.createMassMatrix();
            obj.cost            = Cost(s);
        end

        function M = createMassMatrix(obj)
            test   = LagrangianFunction.create(obj.mesh, 1, 'P1');
            trial  = LagrangianFunction.create(obj.mesh, 1, 'P1');
            f = @(u,v) DP(v,u);
            M = IntegrateLHS(f,test,trial,obj.mesh,'Domain',2);
            M = diag(sum(M,1));
        end

        function createConstraint(obj)
            s.shapeFunctions{1} = obj.volume;
            s.Msmooth           = obj.createMassMatrix();
            obj.constraint      = Constraint(s);
        end

       function createDualVariable(obj)
            s.nConstraints   = 1;
            l                = DualVariable(s);
            obj.dualVariable = l;
        end

        function createOptimizer(obj)
            s.monitoring     = true;
            s.cost           = obj.cost;
            s.constraint     = obj.constraint;
            s.designVariable = obj.designVariable;
            s.dualVariable   = obj.dualVariable;
            s.maxIter        = 1000;
            s.tolerance      = 1e-8;
            s.constraintCase = {'EQUALITY'};
            s.ub             = 1;
            s.lb             = 0;
            s.volumeTarget   = 0.4;
            s.primal         = 'PROJECTED GRADIENT';
            s.gif            = true;
            s.gifName        = 'Tutorial05_1';
            s.printing       = true;
            s.printName      = 'Tutorial05_1';
            opt              = OptimizerMMA(s);
            opt.solveProblem();
            obj.optimizer = opt;
        end

        function m = createMaterial(obj)
            x = obj.designVariable;
            f = x.obtainDomainFunction();
            f = f{1}.project('P1'); 
            s.density  =  f;
            s.type     = 'HomogenizedMicrostructure';
            s.mesh     = obj.mesh;
            s.young    = 1.0;
            s.fileName = 'HomogenizationResultsReinforcedHexagon';
            m = MaterialFactory.create(s);

        end

          function bc = createBoundaryConditions(obj)
            xMin = min(obj.mesh.coord(:,1));
            xMax = max(obj.mesh.coord(:,1));
            yMin = min(obj.mesh.coord(:,2));
        
            
            isBottom = @(coor) abs(coor(:,2) - yMin) < 1e-12;
        
            
            isDirLeft  = @(coor) isBottom(coor) & ...
                                 coor(:,1) >= xMin        & coor(:,1) <= xMin + 0.2;
            isDirRight = @(coor) isBottom(coor) & ...
                                 coor(:,1) >= xMax - 0.2 & coor(:,1) <= xMax;
        
            
            isForce = @(coor) isBottom(coor) & ...
                              coor(:,1) >= xMin + 0.9 & coor(:,1) <= xMin + 1.1;
        
            
            sDir{1}.domain    = @(coor) isDirLeft(coor);
            sDir{1}.direction = [1,2];
            sDir{1}.value     = 0;
        
            sDir{2}.domain    = @(coor) isDirRight(coor);
            sDir{2}.direction = [1,2];
            sDir{2}.value     = 0;
        
            
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