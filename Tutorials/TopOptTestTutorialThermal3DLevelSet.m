classdef TopOptTestTutorialThermal3DLevelSet < handle

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
        thermalCompliance
        source
        primalUpdater
    end

    methods (Access = public)

        function obj = TopOptTestTutorialThermal3DLevelSet()
            for f0 = [1e-2, 1e-3, 1e-4]
                obj.init()
                obj.createMesh();
                obj.createDesignVariable();
                obj.createFilter();
                obj.createMaterialInterpolator(f0);
                obj.createThermalProblem();
                obj.createThermalCompliance();
                obj.createVolumeConstraint();
                obj.createCost();
                obj.createConstraint();
                obj.createDualVariable();
                obj.createPrimalUpdater();
                obj.createOptimizer(f0);
            end
        end

    end

    methods (Access = private)

        function init(obj)
            close all;
        end

        function createMesh(obj)
            obj.mesh = HexaMesh(1,1,1,30,30,30);
        end

        function createDesignVariable(obj)

            s.type = 'Full';
            g      = GeometricalFunction(s);
            lsFun  = g.computeLevelSetFunction(obj.mesh);
            s.fun  = lsFun;
            s.mesh = obj.mesh;
            s.type = 'LevelSet';
            s.plotting = true;
            ls     = DesignVariable.create(s);           
            obj.designVariable = ls;
        end

        function createFilter(obj)
            s.filterType = 'PDE';
            s.mesh  = obj.mesh;
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            f = Filter.create(s);
            f.updateEpsilon(2*obj.mesh.computeMinCellSize())
            obj.filter = f;
        end

        function createMaterialInterpolator(obj,f0) % Conductivity
%             s.interpolation  = 'SIMPThermal';   
            s.interpolation  = 'SimpAllThermal';
            s.f0   = f0;                                             
            s.f1   = 1;  
            s.dim  = '3D'
%             s.pExp = 3;
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
            s.dim = '3D';
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

        function uMesh = createBaseDomain(obj)
            sG.type          = 'Full';
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
            s.uMesh  = obj.createBaseDomain();
            s.filter = obj.filter;
            s.test = LagrangianFunction.create(obj.mesh,1,'P1');
            s.volumeTarget = 0.2;
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
            n = obj.mesh.nnodes;
            h = obj.mesh.computeMinCellSize();
            M = h^2*sparse(1:n,1:n,ones(1,n),n,n);
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

        function createPrimalUpdater(obj)
            s.mesh = obj.mesh;
            obj.primalUpdater = SLERP(s);
        end

         function createOptimizer(obj,f0)
            s.shallPrint     = true;
            s.monitoring     = true;
            s.cost           = obj.cost;
            s.constraint     = obj.constraint;
            s.designVariable = obj.designVariable;
            s.primalUpdater     = obj.primalUpdater;
            s.GIFname        = string(f0)+'-';
            s.maxIter        = 1000;
            s.tolerance      = 1e-8;
            s.constraintCase{1} = 'EQUALITY';
            s.primal         = 'SLERP';
            s.ub             = inf;
            s.lb             = -inf;
            s.etaNorm        = 0.02; % 0.5
            s.etaNormMin     = 0.02;
            s.gJFlowRatio    = 1.0; %5.0;    %0.2  2.0; 60.0
            s.etaMax         = 40.0; %0.1;    % 1 - 5.0 5.0
            s.etaMaxMin      = 0.02; %0.05; 
            s.filter         = obj.filter;
            opt = OptimizerNullSpace(s);
            opt.solveProblem();
            obj.optimizer = opt;
            saveas(figure(1),string(f0)+'design.png','png')
            saveas(figure(2),string(f0)+'graficos.png','png')
            writematrix(obj.designVariable.fun.fValues,string(f0)+'.txt')
            obj.designVariable.fun.print('level_set'+string(f0))
        end


        function bc = createBoundaryConditions(obj)
            xMax    = max(obj.mesh.coord(:,1));
            yMax    = max(obj.mesh.coord(:,2));
            zMin    = min(obj.mesh.coord(:,3));
            isDir   = @(coor) abs(coor(:,3))==zMin & abs(coor(:,1))>=0.45*xMax & abs(coor(:,1))<=0.55*xMax & abs(coor(:,2))>=0.45*yMax & abs(coor(:,2))<=0.55*yMax;  
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