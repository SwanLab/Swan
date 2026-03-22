classdef TutorialCohesive < handle

    properties (Access = public)
        output  
    end

    properties (Access = private)
        mesh
        boundaryConditions
        functional
        solverType
        cohesiveMesh
        tractionSeparation
    end

    methods (Access = public)

        function obj = TutorialCohesive()
            obj.init();
            obj.createMesh();
            obj.defineCase(); %Benchmarks
            obj.createCohesiveFunctional()
            obj.solveProblem()
        end

        function solveProblem(obj)
            s.mesh = obj.cohesiveMesh.mesh;
            s.boundaryConditions = obj.boundaryConditions;
            s.functional = obj.functional;

            s.tolerance              = 1e-8;
            s.maxIter                = 20;

            s.solverType             = obj.solverType;
            
            CohComp = CohesiveComputer(s);
            CohComp.compute();
            obj.output = CohComp.data;
        end

    end

    methods (Access = private)

        function init(obj)
            close all
        end

        function createMesh(obj)
            s.baseMesh = UnitQuadMesh(3,3);
            ymin = 0;
            s.isFractured = @(coord) abs(coord(:,2) - ymin) <=1e-10;
            obj.cohesiveMesh = CohesiveMesh(s);
        end

        function defineCase(obj)
            obj.solverType = 'Newton';
            obj.createBoundaryConditions();
            obj.createTractionSeparation();
        end

        function createBoundaryConditions(obj)
            bc.type = 'ForceTractionY';
            bc.values = [0:0.1:20];
            obj.boundaryConditions  = BoundaryConditionsCreator(obj.cohesiveMesh.mesh,bc);
        end

        function createTractionSeparation(obj)
            s.K         = 1e10;
            s.jumpCrit  = 0.1;
            s.jumpFinal = 0.01;
            s.lawType = 'TractionBiliniarUncoupled';
            obj.tractionSeparation = CohesiveTractionSeparation(s);
        end

        function createCohesiveFunctional(obj)
            s.tractionSeparation = obj.tractionSeparation;
            s.material           = obj.createMaterial();
            s.mesh               = obj.cohesiveMesh.mesh;
            s.cohesiveMesh       = obj.cohesiveMesh;
            s.quadOrder          = 2;
            s.test = LagrangianFunction.create(obj.cohesiveMesh.mesh,2,'P1');
            s.cohTest = LagrangianFunction.create(obj.cohesiveMesh.lineMesh,2,'P1');
            obj.functional = CohesiveFunctional(s);
        end

        function m = createMaterial(obj)
            k.type    = 'ISOTROPIC';
            k.ptype   = 'ELASTIC';
            k.ndim    = obj.cohesiveMesh.mesh.ndim;
            k.young   = ConstantFunction.create(100,obj.cohesiveMesh.mesh);
            k.poisson = ConstantFunction.create(0.33, obj.cohesiveMesh.mesh);
            k.cohesiveMesh = obj.cohesiveMesh;
            m    = Material.create(k);
        end

    end

end