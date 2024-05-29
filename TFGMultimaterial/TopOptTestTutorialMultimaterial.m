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
    end

    methods (Access = public)

        function obj = TopOptTestTutorialMultimaterial()
            obj.init()
            obj.createMesh();
            obj.createDesignVariable();
            obj.createFilter();
            obj.createInterpolators();
            obj.createBoundaryConditions();
            obj.createElasticProblem();
            obj.createComplianceFromConstiutive();
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
            obj.nMat = 4;
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

            s.type                 = 'LevelSet';
            s.plotting             = false;
            s.fValues              = lsFun{1};
            s.order                = 'P1';
            s.fun                  = LagrangianFunction(s);
            ls1                    = DesignVariable.create(s);
            obj.designVariable{1} = ls1;

            s.fValues              = lsFun{2};
            s.fun                  = LagrangianFunction(s);
            ls1                    = DesignVariable.create(s);
            obj.designVariable{2} = ls1;

            s.fValues              = lsFun{3};
            s.fun                  = LagrangianFunction(s);
            ls1                    = DesignVariable.create(s);
            obj.designVariable{3} = ls1;
        end

        function createFilter(obj)
            s.filterType = 'LUMP';
            s.mesh  = obj.mesh;
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            f = Filter.create(s);
            obj.filter = f;
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

         function createVolumeConstraint(obj)
            s.mesh   = obj.mesh;
            s.filter = obj.filter;
            s.gradientTest = LagrangianFunction.create(obj.mesh,1,'P1');
            s.volumeTarget = [0.2 0.2 0.2 1.2]; % prova
            v = VolumeConstraint(s);
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
            s.maxIter        = 100;
            s.tolerance      = 1e-8;
            s.constraintCase = {'EQUALITY'};
            s.volumeTarget   = [0.2 0.2 0.2 1.2]; %prova
            s.primal         = 'SLERP';
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
            
            constituitiveTensor = ElasticTensorComputer(s);
            m = constituitiveTensor.C;
        end

        function M = createMassMatrix(obj)
            s.test  = LagrangianFunction.create(obj.mesh,1,'P1');
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            s.mesh  = obj.mesh;
            s.type  = 'MassMatrix';
            LHS = LHSintegrator.create(s);
            M = LHS.compute;
        end


    end
end