classdef TutorialCohesive < handle

    properties (Access = public)
        output  
    end

    properties (Access = private)
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
            s.mesh = obj.cohesiveMesh.fullMesh;
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
            s.baseMesh = UnitQuadMesh(2,2);
            y = 0;
            s.isFractured = @(coord) abs(coord(:,2) - y) <= 1e-10;
            obj.cohesiveMesh = CohesiveMesh(s);
        end

        function defineCase(obj)
            obj.solverType = 'Newton';
            obj.createBoundaryConditions();
            obj.createTractionSeparation();
        end

        function createBoundaryConditions(obj)
            bc.type = 'DisplacementTractionY';
            bc.values = 0.05;
            obj.boundaryConditions  = BoundaryConditionsCreator(obj.cohesiveMesh.fullMesh,bc);
        end

        function createTractionSeparation(obj)
            s.jumpCrit  = 0.001;
            s.jumpFinal = 0.2;
            s.fractureStrength  = 1;
            s.fractureToughness = 1;
            s.lawType = 'TractionBiliniarCoupled';
            % s.lawType = 'TractionBiliniarUncoupled';
            obj.tractionSeparation = CohesiveTractionSeparation(s);
        end

        function createCohesiveFunctional(obj)
            s.tractionSeparation = obj.tractionSeparation;
            s.material           = obj.createMaterial();
            s.mesh               = obj.cohesiveMesh.fullMesh;
            s.cohesiveMesh       = obj.cohesiveMesh;
            s.quadOrder          = 2;
            s.test = LagrangianFunction.create(obj.cohesiveMesh.fullMesh,2,'P1');
            obj.functional = CohesiveFunctional(s);
        end

        function m = createMaterial(obj)
            k.type    = 'ISOTROPIC';
            k.ptype   = 'ELASTIC';
            k.ndim    = obj.cohesiveMesh.fullMesh.ndim;
            k.young   = ConstantFunction.create(1e6,obj.cohesiveMesh.fullMesh);
            k.poisson = ConstantFunction.create(0, obj.cohesiveMesh.fullMesh);
            k.cohesiveMesh = obj.cohesiveMesh;
            m    = Material.create(k);
        end

    end

end