classdef IterativeTutorialCohesive < handle

    properties (Access = public)
        output
        inputData
    end

    properties (Access = private)
        boundaryConditions
        functional
        solverType
        cohesiveMesh
        tractionSeparation
    end

    methods (Access = public)

        function obj = IterativeTutorialCohesive(inputData)

            obj.inputData = inputData;

            obj.createMesh();
            obj.defineCase();
            obj.createCohesiveFunctional();
            obj.solveCohesiveProblem();
        end

        function solveCohesiveProblem(obj)

            s.cohesiveMesh       = obj.cohesiveMesh;
            s.boundaryConditions = obj.boundaryConditions;
            s.functional         = obj.functional;
            s.tolerance          = 1e-4;
            s.maxIter            = 500;
            s.tractionLaw        = obj.tractionSeparation;

            s.monitoring.set     = true;
            s.monitoring.print   = true;

            s.solverType         = obj.solverType;

            CohComp = CohesiveComputer(s);
            CohComp.compute();

            obj.output = CohComp.data;

        end

    end

    methods (Access = private)

        function createMesh(obj)

            s.xCohLineMax = obj.inputData.xCohLineMax;
            s.yCohLine    = obj.inputData.yCohLine;

            s.baseMesh = QuadMesh( ...
                obj.inputData.l,...
                obj.inputData.h,...
                obj.inputData.nx,...
                obj.inputData.ny);

            obj.cohesiveMesh = NewCohesiveMesh(s);

        end

        function defineCase(obj)

            obj.solverType = obj.inputData.solverType;

            obj.createBoundaryConditions();
            obj.createTractionSeparation();

        end

        function createBoundaryConditions(obj)

            bc.type   = obj.inputData.problemType;
            bc.values = obj.inputData.bcValues;

            obj.boundaryConditions = BoundaryConditionsCreator( obj.cohesiveMesh.fullMesh,bc);

        end

        function createTractionSeparation(obj)

            obj.tractionSeparation = ...
                CohesiveTractionSeparation(obj.inputData);

        end

        function createCohesiveFunctional(obj)

            s.tractionSeparation = obj.tractionSeparation;
            s.material           = obj.createMaterial();
            s.mesh               = obj.cohesiveMesh.fullMesh;
            s.cohesiveMesh       = obj.cohesiveMesh;
            s.quadOrder          = 2;
            s.test = LagrangianFunction.create(obj.cohesiveMesh.fullMesh, 2, 'P1');

            obj.functional = CohesiveFunctional(s);

        end

        function m = createMaterial(obj)

            k.type    = 'ISOTROPIC';
            k.ptype   = 'ELASTIC';
            k.ndim    = obj.cohesiveMesh.fullMesh.ndim;

            k.young = ConstantFunction.create( ...
                obj.inputData.young,...
                obj.cohesiveMesh.fullMesh);

            k.poisson = ConstantFunction.create( ...
                obj.inputData.poisson,...
                obj.cohesiveMesh.fullMesh);

            k.cohesiveMesh = obj.cohesiveMesh;

            m = Material.create(k);
            m.setCohesiveMesh(obj.cohesiveMesh);

        end

    end

end