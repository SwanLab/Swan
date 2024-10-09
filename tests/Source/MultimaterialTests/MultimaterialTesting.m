classdef MultimaterialTesting < handle

    properties (Access = private)
        mesh
        designVariable
        area
        filter
        physicalProblem
        compliance
        volume
        cost
        constraint
        dualVariable
        optimizer
        nMat
        pdeCoeff
        matInterp
        mat
        bc
        nLevelSet
    end

    methods (Access = public)

        function obj = MultimaterialTesting()
            obj.init()
            obj.createMesh();
            obj.createDesignVariable();
            obj.createInterpolators();
            obj.createBoundaryConditions();
            obj.createElasticProblem();
            obj.createCompliance();
            obj.createVolumeConstraint();
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
            obj.nMat = 4; % including weakest mat
            obj.nLevelSet = 3;
        end

        function createMesh(obj)
            %obj.mesh = TriangleMesh(6,1,150,25); % Bridge
            obj.mesh = TriangleMesh(2,1,20,10); % Beam
            %obj.mesh = TriangleMesh(2,1,100,50); % Arch
            p = obj.mesh.coord';
            t = obj.mesh.connec';
            obj.area = pdetrg(p,t);
        end

        function createDesignVariable(obj)
            s.mesh                 = obj.mesh;
            s.type                 = 'Full';
            lsFun{1}               = @(x) -ones(size(x(1,:,:)));
            lsFun{2}               = @(x) -cos(x(1,:,:))+0.5;
            lsFun{3}               = @(x) sin(x(1,:,:))-0.5;
            
            s.type                 = 'MultiLevelSet';
            s.lsFun                = lsFun;
            s.mesh                 = obj.mesh;
            s.unitM                = obj.createMassMatrix();
            s.plotting             = true;
            obj.designVariable     = DesignVariable.create(s);
        end

        function createInterpolators(obj)
            matProp      = MaterialPropertiesComputer(); 
            obj.mat.A    = matProp.matA; % MATERIAL 1
            obj.mat.B    = matProp.matB; % MATERIAL 2
            obj.mat.C    = matProp.matC; % MATERIAL 3
            obj.mat.D    = matProp.matD; % VOID
            s.mat        = obj.mat;
            s.m          = obj.mesh;
            
            obj.pdeCoeff = PDECoefficientsComputer(s);  


            sC.E = [obj.mat.A.young,obj.mat.B.young,obj.mat.C.young,obj.mat.D.young];
            sC.C0 = obj.pdeCoeff.tensor(:,1);
            sC.nu1 = obj.mat.A.nu;
            obj.matInterp = MultiMaterialInterpolation(sC);
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
            obj.bc = BoundaryConditions(s);
        end

        function createElasticProblem(obj)
            s.mesh = obj.mesh;
            s.scale = 'MACRO';
            s.material = obj.createMaterial();
            s.dim = '2D';
            s.boundaryConditions = obj.bc;
            %s.interpolationType = 'LINEAR';
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
            s.nMat = obj.nMat;
            s.mat = obj.mat;
            s.mesh = obj.mesh;
            s.pdeCoeff = obj.pdeCoeff;
            s.bc = obj.bc;
            s.complianceFromConstitutive = obj.createComplianceFromConstiutive();
            s.material = obj.createMaterial();
            s.materialInterpolator = obj.matInterp;
            c = ComplianceFunctionalComputer(s);
            obj.compliance = c;
        end

         function createVolumeConstraint(obj)
            s.volumeTarget = [0.2 0.2 0.2 1.4]; 
            s.mesh = obj.mesh;
            s.area  = obj.area;
            v = VolumeConstraintComputer(s);
            obj.volume = v;
         end

         function createCost(obj)
            s.shapeFunctions{1} = obj.compliance;
            s.weights           = 1;
            s.Msmooth           = obj.createMassMatrix();
            obj.cost            = Cost(s);
         end

         function createConstraint(obj)
            s.shapeFunctions{1} = obj.volume;
            s.Msmooth           = obj.createMassMatrix();
            obj.constraint      = Constraint(s);
         end

         function createDualVariable(obj)
            s.nConstraints   = obj.nMat-1;
            l                = DualVariable(s);
            obj.dualVariable = l;
         end

         function createOptimizer(obj)
            s.monitoring     = true;
            s.cost           = obj.cost;
            s.constraint     = obj.constraint;
            s.designVariable = obj.designVariable;
            s.dualVariable   = obj.dualVariable;
            s.maxIter        = 2;
            s.tolerance      = 1e-8;
            s.constraintCase = repmat({'EQUALITY'},[obj.nMat,1]);
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
            s.materialInterpolator = obj.matInterp;
            s.dim                  = '2D';
            m = Material.create(s);
        end

        function M = createMassMatrix(obj)
            nnodes  = obj.mesh.nnodes*obj.nLevelSet;
            %nnodes  = obj.mesh.nnodes;
            indices = transpose(1:nnodes);
            vals    = ones(size(indices));
            h       = obj.mesh.computeMeanCellSize();
            M       = h^2*sparse(indices,indices,vals,nnodes,nnodes);
        end


    end
end