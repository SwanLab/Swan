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
            
            problem      = 'DisplacementTractionY';
            nameMesh     = 'M10by4';
            obj.init();
            obj.createMesh(problem,nameMesh);
            obj.defineCase(problem);
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
    
        function createMesh(obj,problem,nameMesh)
            s = obj.makeBenchmarkMesh(problem,nameMesh);
            s.baseMesh       = UnitQuadMesh(1,1);
            obj.cohesiveMesh = CohesiveMesh(s);
        end

        function defineCase(obj,problem)
            obj.solverType = 'Newton';
            obj.createBoundaryConditions(problem);
            obj.createTractionSeparation();
        end

        function createBoundaryConditions(obj,problem)
            bc.type = problem;
            % bc.type = 'DisplacementTractionY';
            bc.values = [0.15];
            obj.boundaryConditions  = BoundaryConditionsCreator(obj.cohesiveMesh.fullMesh,bc);
        end

        function createTractionSeparation(obj)
            s.jumpCrit  = 0.001;
            s.jumpFinal = 0.2;
            s.fractureStrength  = 1;
            s.fractureToughness = 1;
            % s.lawType = 'TractionBiliniarCoupled';
            s.lawType = 'TractionBiliniarUncoupled';
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
            k.poisson = ConstantFunction.create(0.3, obj.cohesiveMesh.fullMesh);
            k.cohesiveMesh = obj.cohesiveMesh;
            m    = Material.create(k);
            m.setCohesiveMesh(obj.cohesiveMesh);
        end

        function c = makeBenchmarkMesh(obj, problem,nameMesh)
            a.fileName = nameMesh;
            s = FemDataContainer(a);
            c.baseMesh = s.mesh;

            switch problem
                case 'DoubleCantileverBeam'
                    xmax = 72*0.001; y = 3.12 * 0.5 * 0.001;
                case 'DisplacementTractionY'
                    
            end

            c.isFracturedLine  = @(coord) abs(coord(:,2) - y) <= 1e-10;
            c.isFracturedUntil = @(coord) coord(:,1) - xmax   <= 1e-10;
        end

    end

end