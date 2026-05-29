classdef TutorialDehomoFirst < handle

    properties (Access = private)
        mesh
        filterCompliance
        %designVariable
        materialInterpolator
        physicalProblem
        compliance
        cost
        optimizer
        filterRegularization
    end
    properties (Access = public)
        designVariable
    end

    methods (Access = public)

        function obj = TutorialDehomoFirst()
            obj.init()
            obj.createMesh();
            obj.createDesignVariable();
            obj.createFilterRegularization();
            obj.createElasticProblem();
            obj.createBaseDomain();
            obj.createComplianceFromConstiutive();
            obj.createCost();
            obj.createOptimizer();
        end
        
       
    end

    methods (Access = private)

        function init(obj)
            close all;
        end

        function createMesh(obj)
            
            obj.mesh = TriangleMesh(2,1,100,100);
        end

        function createDesignVariable(obj)
            s.fHandle = @(x) zeros(size(x(1,:,:)));
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

       
        function createFilterRegularization(obj)
            eOverhmin = 6;
            s.filterType    = 'PDE';
            s.mesh          = obj.mesh;
            s.boundaryType  = 'Neumann';  
            s.metric        = 'Isotropy';
            s.trial         = LagrangianFunction.create(obj.mesh, 1, 'P1');
            f = Filter.create(s);
            
            epsilon = eOverhmin * obj.mesh.computeMeanCellSize();
            f.updateEpsilon(epsilon);
            
            obj.filterRegularization = f;
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
        function c = createCompliance(obj)
            s.mesh                        = obj.mesh;
            s.filter                      = obj.filterRegularization;
            s.complainceFromConstitutive  = obj.createComplianceFromConstiutive();
            s.material                    = obj.createMaterial();
            c = ComplianceFunctional(s);
            
        end


        
        function createCost(obj)
            s.shapeFunctions{1} = obj.createCompliance();
            s.weights           = 1;
            s.Msmooth           = obj.createMassMatrix();
            obj.cost            = Cost(s);
        end

        function M = createMassMatrix(obj)
            test  = LagrangianFunction.create(obj.mesh,1,'P1');
            trial = LagrangianFunction.create(obj.mesh,1,'P1'); 
            M = IntegrateLHS(@(u,v) DP(v,u),test,trial,obj.mesh,'Domain');
        end


        function createOptimizer(obj)
            s.cost           = obj.cost;
            s.designVariable = obj.designVariable;
            s.monitoring     = true;
            s.lb             = -0.6;
            s.ub             = 0.6;
            s.maxIter        = 200;
            opt              = OptimizerProjectedGradient(s);
            opt.solveProblem();
            obj.optimizer = opt;
            

            

            
        end

    
        function m = createMaterial(obj)
            x = obj.designVariable;
            
            
            s.density  = x;
            s.type     = 'HomogenizedMicrostructure';
            s.mesh     = obj.mesh;
            s.young    = 1.0;
            s.fileName = 'HomogenizationLattice4';
            m = MaterialFactory.create(s);
        end

        function bc = createBoundaryConditions(obj)
            
            xMin = min(obj.mesh.coord(:,1));
            xMax = max(obj.mesh.coord(:,1));
            yMin = min(obj.mesh.coord(:,2));
            yMax = max(obj.mesh.coord(:,2));
            
            
            isBottom = @(coor) abs(coor(:,2) - yMin) < 1e-12;
            
            
            isDirLeft = @(coor) isBottom(coor) & ...
                                coor(:,1) >= xMin & coor(:,1) <= xMin + 0.2*xMax;
            
            
            isDirRight = @(coor) isBottom(coor) & ...
                                 coor(:,1) >= xMax - 0.2*xMax & coor(:,1) <= xMax;
            
            
            isForce = @(coor) isBottom(coor) & ...
                               coor(:,1) >= 0.45*xMax & coor(:,1) <= 0.55*xMax;
            
            
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
            
            s.periodicFun = [];
            s.mesh = obj.mesh;
            bc = BoundaryConditions(s);
        end
    end
end