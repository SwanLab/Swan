classdef MultimaterialTesting < handle

    properties (Access = private)
        mesh
        designVariable
        materialInterpolator
        filter
        boundaryConditions
        physicalProblem
        compliance
        volumeA
        volumeB
        volumeC
        cost
        constraint
        dualVariable
        optimizer
    end

    methods (Access = public)

        function obj = MultimaterialTesting()
            obj.init()
            obj.createMesh();
            obj.createFilter();
            obj.createDesignVariable();
            obj.createMaterialInterpolator();
            obj.createBoundaryConditions();
            obj.createElasticProblem();
            obj.createComplianceFromConstiutive();
            obj.createCompliance();
            obj.createVolumeConstraints();
            obj.createCost();
            obj.createConstraint();
            obj.createDualVariable();
            obj.createOptimizer();
        end

        function lsVals = solve(obj)
            obj.optimizer.solveProblem();
            lsVals = obj.designVariable.fun.fValues;
        end

    end

    methods (Access = private)

        function init(obj)
            close all;
        end

        function createMesh(obj)
            obj.mesh = TriangleMesh(2,1,14,7);
        end

        function createFilter(obj)
            s.filterType = 'LUMP';
            s.mesh       = obj.mesh;
            s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
            f            = Filter.create(s);
            obj.filter   = f;
        end

        function createDesignVariable(obj)
            lsFun{1} = obj.createLevelSetFunction(@(x) -ones(size(x(1,:,:))));
            lsFun{2} = obj.createLevelSetFunction(@(x) -cos(x(1,:,:))+0.5);
            lsFun{3} = obj.createLevelSetFunction(@(x) sin(x(1,:,:))-0.5);
            
            s.type             = 'MultiLevelSet';
            s.lsFun            = lsFun;
            s.mesh             = obj.mesh;
            s.plotting         = true;
            obj.designVariable = DesignVariable.create(s);
        end

        function lsFun = createLevelSetFunction(obj,fH)
            s.type    = 'Given';
            s.fHandle = fH;
            g         = GeometricalFunction(s);
            lsFun     = g.computeLevelSetFunction(obj.mesh);
        end
        
        function createMaterialInterpolator(obj)
            s.E    = [1,0.5,0.25,1e-3];
            s.nu   = (1/3)*[1,1,1,1];
            s.ndim = 2;
            obj.materialInterpolator = MultiMaterialInterpolation(s);
        end

        function createBoundaryConditions(obj)
            xMax    = max(obj.mesh.coord(:,1));
            yMax    = max(obj.mesh.coord(:,2));
            isDir   = @(coor)  abs(coor(:,1))==0;
            isForce = @(coor)  (abs(coor(:,1))==xMax & abs(coor(:,2))>=0.35*yMax & abs(coor(:,2))<=0.65*yMax);

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
            for i = 1:numel(sPL)
                pl = PointLoad(obj.mesh, sPL{i});
                pointloadFun = [pointloadFun, pl];
            end
            s.pointloadFun = pointloadFun;

            s.periodicFun  = [];
            s.mesh = obj.mesh;
            obj.boundaryConditions = BoundaryConditions(s);
        end

        function createElasticProblem(obj)
            s.mesh = obj.mesh;
            s.scale = 'MACRO';
            s.material = obj.createMaterial();
            s.dim = '2D';
            s.boundaryConditions = obj.boundaryConditions;
            s.solverType = 'REDUCED';
            s.solverMode = 'DISP';
            s.solverCase = 'DIRECT';
            fem = ElasticProblem(s);
            obj.physicalProblem = fem;
        end

        function c = createComplianceFromConstiutive(obj)
            s.mesh         = obj.mesh;
            s.stateProblem = obj.physicalProblem;
            c = ComplianceFromConstiutiveTensor(s);
        end

        function createCompliance(obj)
            s.mesh                       = obj.mesh;
            s.filter                     = obj.filter;
            s.complainceFromConstitutive = obj.createComplianceFromConstiutive();
            s.material                   = obj.createMaterial();
            c = ComplianceFunctional(s);
            obj.compliance = c;
        end

        function createVolumeConstraints(obj)
            obj.volumeA = obj.createIndivVolumeConstraint(0.1,1);
            obj.volumeB = obj.createIndivVolumeConstraint(0.1,2);
            obj.volumeC = obj.createIndivVolumeConstraint(0.1,3);
        end

        function v = createIndivVolumeConstraint(obj,target,ID)
            s.volumeTarget = target;
            s.nMat         = 4;
            s.matID        = ID;
            s.mesh         = obj.mesh;
            s.gradientTest = LagrangianFunction.create(obj.mesh,1,'P1');
            v              = MultiMaterialVolumeConstraint(s);
         end

         function createCost(obj)
            s.shapeFunctions{1} = obj.compliance;
            s.weights           = 1;
            s.Msmooth           = obj.createMassMatrix();
            obj.cost            = Cost(s);
         end

         function createConstraint(obj)
            s.shapeFunctions{1} = obj.volumeA;
            s.shapeFunctions{2} = obj.volumeB;
            s.shapeFunctions{3} = obj.volumeC;
            s.Msmooth           = obj.createMassMatrix();
            obj.constraint      = Constraint(s);
         end

         function createDualVariable(obj)
            s.nConstraints   = 3;
            l                = DualVariable(s);
            obj.dualVariable = l;
         end

         function createOptimizer(obj)
            s.monitoring     = false;
            s.cost           = obj.cost;
            s.constraint     = obj.constraint;
            s.designVariable = obj.designVariable;
            s.dualVariable   = obj.dualVariable;
            s.maxIter        = 1;
            s.tolerance      = 1e-8;
            s.constraintCase = repmat({'EQUALITY'},[3,1]);
            s.primal         = 'SLERP';
            s.ub             = inf;
            s.lb             = -inf;
            s.etaNorm        = inf;
            s.gJFlowRatio    = 2;
            obj.optimizer    = OptimizerNullSpace(s);
        end

        function m = createMaterial(obj)
            x                      = obj.designVariable;
            s.type                 = 'DensityBased';
            s.density              = x;
            s.materialInterpolator = obj.materialInterpolator;
            s.dim                  = '2D';
            s.mesh                 = obj.mesh;
            m = Material.create(s);
        end

        function M = createMassMatrix(obj)
            nnodes  = obj.mesh.nnodes*3;
            indices = transpose(1:nnodes);
            vals    = ones(size(indices));
            h       = obj.mesh.computeMeanCellSize();
            M       = h^2*sparse(indices,indices,vals,nnodes,nnodes);
        end


    end
end