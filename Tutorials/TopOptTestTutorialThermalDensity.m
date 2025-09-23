classdef TopOptTestTutorialThermalDensity < handle

    properties (Access = private)
        mesh
        filter
        designVariable
        materialInterpolator
        physicalProblem
        thermalCompliance
        volume
        cost
        constraint
        dualVariable
        optimizer
        source 
    end

    methods (Access = public)

        function obj = TopOptTestTutorialThermalDensity()
            obj.init()
            obj.createMesh();
            obj.createDesignVariable();
            obj.createFilter();
            obj.createMaterialInterpolator();
            obj.createThermalProblem();

            obj.createThermalCompliance();
            obj.createVolumeConstraint();
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
            %UnitMesh better
            x1      = linspace(0,1,50);
            x2      = linspace(0,1,50);
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
            s.fun     = aFun.project('P1');
            s.mesh    = obj.mesh;
            s.type = 'Density';
            s.plotting = true;
            dens    = DesignVariable.create(s);
            obj.designVariable = dens;
        end

        function createFilter(obj)
            s.filterType = 'LUMP';
            s.mesh  = obj.mesh;
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            f = Filter.create(s);
            obj.filter = f;
        end

        function createMaterialInterpolator(obj) % Conductivity
            s.interpolation  = 'SIMPThermal';   
            s.f0   = 0.01;                                             
            s.f1   = 1;                                                    
            s.pExp = 3;
            a = MaterialInterpolator.create(s);
            obj.materialInterpolator = a;            
        end   

        function createThermalProblem(obj)
            s.mesh = obj.mesh;
            s.conductivity = obj.materialInterpolator; 
            Q = LagrangianFunction.create(obj.mesh,1,'P1');
            fValues = ones(Q.nDofs,1);
            Q.setFValues(fValues);
            s.source       = Q;  
            s.dim = '2D';
            s.boundaryConditions = obj.createBoundaryConditions();
            s.interpolationType = 'LINEAR';
            s.solverType = 'REDUCED';
            s.solverMode = 'DISP';
            s.solverCase = 'DIRECT';
            fem = ThermalProblem(s); 
            obj.physicalProblem = fem;
        end

        function createThermalCompliance(obj)
            s.mesh                        = obj.mesh;
            s.filter                      = obj.filter;
            s.stateProblem                = obj.physicalProblem;
            s.conductivity                =  obj.materialInterpolator; 
            c = ThermalComplianceFunctional(s);  
            obj.thermalCompliance = c;
        end

        function createVolumeConstraint(obj)
            s.mesh   = obj.mesh;
            s.filter = obj.filter;
            s.gradientTest = LagrangianFunction.create(obj.mesh,1,'P1');
            s.volumeTarget = 0.4;
            v = VolumeConstraint(s);
            obj.volume = v;
        end

        function createCost(obj)
            s.shapeFunctions{1} = obj.thermalCompliance;
            s.weights           = 1;
            s.Msmooth           = obj.createMassMatrix();
            obj.cost            = Cost(s);
        end

        function M = createMassMatrix(obj)
            s.test  = LagrangianFunction.create(obj.mesh,1,'P1');
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            s.mesh  = obj.mesh;
            s.type  = 'MassMatrix';
            LHS = LHSIntegrator.create(s);
            M = LHS.compute;     
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
            s.maxIter        = 700;
            s.tolerance      = 1e-8;
            s.constraintCase = {'EQUALITY'};
            s.ub             = 1;
            s.lb             = 0;
            s.volumeTarget   = 0.4;
            s.primal         = 'PROJECTED GRADIENT';
            opt              = OptimizerMMA(s);
            opt.solveProblem();
            obj.optimizer = opt;
        end

        function bc = createBoundaryConditions(obj)
            yMin    = min(obj.mesh.coord(:,2));
            xMax    = max(obj.mesh.coord(:,1));
            isDir   = @(coor) abs(coor(:,2))==yMin & abs(coor(:,1))>=0.4*xMax & abs(coor(:,1))<=0.6*xMax;  
            
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
            bc = BoundaryConditions(s);  
        end
    end
end
