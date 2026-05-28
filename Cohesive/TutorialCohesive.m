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
            obj.solveCohesiveProblem()
        end

        function solveCohesiveProblem(obj)
            s.cohesiveMesh           = obj.cohesiveMesh;
            s.boundaryConditions     = obj.boundaryConditions;
            s.functional             = obj.functional;
            s.tolerance              = 1e-8;
            s.maxIter                = 50;
            s.tractionLaw            = obj.tractionSeparation;
            
            s.monitoring.set         = true;
            s.monitoring.print       = true;

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
            obj.cohesiveMesh = NewCohesiveMesh(s);
        end

        function defineCase(obj,problem)
            obj.solverType = 'Newton';
            obj.createBoundaryConditions(problem);
            obj.createTractionSeparation();
        end

        function createBoundaryConditions(obj,problem)
            bc.type = problem;
            bc.values = 0:0.001:0.5;
            obj.boundaryConditions  = BoundaryConditionsCreator(obj.cohesiveMesh.fullMesh,bc);
        end

        function createTractionSeparation(obj)
            s.tau0Normal          = 15e6;
            s.tau0Shear           = 30e6;
            s.firstCritEnergy     = 260;
            s.secondCritEnergy    = 1002;
            s.eta                 = 2;
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
            k.poisson = ConstantFunction.create(0.3, obj.cohesiveMesh.fullMesh);
            k.cohesiveMesh = obj.cohesiveMesh;
            m    = Material.create(k);
            m.setCohesiveMesh(obj.cohesiveMesh);
        end

        function c = makeBenchmarkMesh(obj, problem,nameMesh)
            switch problem
                case 'DoubleCantileverBeam'
                    c.xCohLineMax = 0.0072; c.yCohLine = 3.12 * 0.5 * 0.001;
                    a.fileName = nameMesh;
                    s = FemDataContainer(a);
                    c.baseMesh = s.mesh;
                case 'DisplacementTractionY' 
                    c.xCohLineMax = 10; c.yCohLine = 0; c.baseMesh = UnitQuadMesh(1,1);
                case 'DisplacementMixed'
                    c.xCohLineMax = 10000; c.yCohLine = 0; c.baseMesh = UnitQuadMesh(1,1);
            end

        end

    end

end