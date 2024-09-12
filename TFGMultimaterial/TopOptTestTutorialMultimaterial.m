classdef TopOptTestTutorialMultimaterial < handle

    properties (Access = private)
        mesh
        area
        filter
        designVariable
        physicalProblem
        compliance
        volume
        cost
        constraint
        dualVariable
        optimizer
        nMat
        pdeCoeff
        mat
        bc
        energy0
        nLevelSet
    end

    methods (Access = public)

        function obj = TopOptTestTutorialMultimaterial()
            obj.init()
            obj.createMesh();
            obj.createDesignVariable();
            obj.createInterpolators();
            obj.createBoundaryConditions();
            obj.createElasticProblem();
            obj.computeEnergy0();
            obj.createCompliance();
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
            obj.nMat = 4; % including weakest mat
            obj.nLevelSet = 3;
        end

        function createMesh(obj)
            %obj.mesh = TriangleMesh(6,1,150,25); % Bridge
            obj.mesh = TriangleMesh(2,1,100,50); % Beam
            %obj.mesh = TriangleMesh(2,1,100,50); % Arch
            p = obj.mesh.coord';
            t = obj.mesh.connec';
            obj.area = pdetrg(p,t);
        end

        function createDesignVariable(obj)
            s.mesh                 = obj.mesh;
            s.type                 = 'Full';
            lsFun{1}               = -ones(size(s.mesh.coord,1),1);
            lsFun{2}               = ones(size(s.mesh.coord,1),1);
            lsFun{3}               = ones(size(s.mesh.coord,1),1);
            
            s.type                 = 'MultiLevelSet';
            s.plotting             = true;
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
        end

        function createBoundaryConditions(obj)
            s.mesh = obj.mesh;
            BoundCond  = BoundaryConditionsSwan(s);
       
            obj.bc = BoundCond.createBoundaryConditionsTutorialBeam();
            %obj.bcSwan = BoundCond.createBoundaryConditionsTutorialBridge();
            %obj.bcSwan = BoundCond.createBoundaryConditionsTutorialArch();
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

        function computeEnergy0(obj)
            s.C = obj.createMaterial();
            s.mesh = obj.mesh;
            s.bc = obj.bc;  
            energy = ComputeInitialEnergy(s);
            obj.energy0 = energy.e0;
        end


        function createCompliance(obj)
            s.energy0 = obj.energy0; % s ha de resoldre un initial elastic problem
            s.nMat = obj.nMat;
            s.mat = obj.mat;
            s.mesh = obj.mesh;
            s.pdeCoeff = obj.pdeCoeff;
            s.bc = obj.bc;
            c = ComplianceFunctionalComputer(s);
            obj.compliance = c;
            %obj.compliance = c.computeFunctionAndGradient(obj.designVariable);
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
            s.maxIter        = 100;
           % s.volumeTarget = [0.4/3 0.4/3 0.4/3 1.6]; % afegit
            s.tolerance      = 1e-8;
            s.constraintCase = repmat({'EQUALITY'},[obj.nMat,1]);
            s.primal         = 'SLERP';
            s.ub             = inf;
            s.lb             = -inf;
            s.etaNorm        = 1;
            s.gJFlowRatio    = 0.5;
            opt = OptimizerNullSpace(s);
            opt.solveProblem();
            obj.optimizer = opt;
        end

        function m = createMaterial(obj)
            s.matProp           = obj.mat;
            s.pdeCoeff          = obj.pdeCoeff;
            s.bc                = obj.bc;
            s.m                 = obj.mesh;
            s.designVariable    = obj.designVariable;
            s.mesh              = obj.mesh;
            
            constituitiveTensor = ElasticTensorComputer(s);
            m = constituitiveTensor.C;
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