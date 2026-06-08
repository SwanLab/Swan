classdef TutorialCohesive < handle

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

        function obj = TutorialCohesive()
            obj.init();
            obj.createMesh();
            obj.defineCase();
            obj.createCohesiveFunctional()
            obj.solveCohesiveProblem()
        end

        function solveCohesiveProblem(obj)
            s.cohesiveMesh           = obj.cohesiveMesh;
            s.boundaryConditions     = obj.boundaryConditions;
            s.functional             = obj.functional;
            s.tolerance              = 1e-6;
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
           
            obj.inputData.young     = 120e9;
            obj.inputData.poisson   = 0.2;
    
            obj.inputData.lawType = 'TractionBiliniarCoupled';
            
            % TURON
            obj.inputData.tau0Normal       = 15e6;
            obj.inputData.tau0Shear        = 30e6;
            obj.inputData.firstCritEnergy  = 260;
            obj.inputData.secondCritEnergy = 1002;
            obj.inputData.eta              = 2;
                        
            % NO TURON
            obj.inputData.jumpCrit         = 1.25e-7;
            obj.inputData.jumpFinal        = 0.025e-3;

            % UnitElem
            obj.inputData.Kcoh        = 3e12;
            obj.inputData.bcValues    = 0:0.000001:0.4;
            obj.inputData.problemType = 'DisplacementTractionY';
            obj.inputData.l = 1;
            obj.inputData.h = 1;
            obj.inputData.yCohLine = 0;
            obj.inputData.xCohLineMax = 1000;
            obj.inputData.nx = 1;
            obj.inputData.ny = 1;

            % % DCB
            % obj.inputData.Kcoh      = 3e12;
            % obj.inputData.bcValues = 0:0.00001:0.008;
            % obj.inputData.problemType = 'DoubleCantileverBeam';
            % obj.inputData.l = 0.103;
            % obj.inputData.h = 3.12e-3;
            % obj.inputData.yCohLine = 0.5*3.12e-3;
            % obj.inputData.xCohLineMax = 72e-3;
            % obj.inputData.nx = 100;
            % obj.inputData.ny = 10;
        end
    
        function createMesh(obj)
            s.xCohLineMax = obj.inputData.xCohLineMax; s.yCohLine = obj.inputData.yCohLine;
            s.baseMesh = QuadMesh(obj.inputData.l ,obj.inputData.h, obj.inputData.nx, obj.inputData.ny);
            obj.cohesiveMesh = NewCohesiveMesh(s);
        end

        function defineCase(obj)
            obj.solverType = 'Newton';
            obj.createBoundaryConditions();
            obj.createTractionSeparation();
        end

        function createBoundaryConditions(obj)
            bc.type = obj.inputData.problemType;
            bc.values = obj.inputData.bcValues;
            obj.boundaryConditions  = BoundaryConditionsCreator(obj.cohesiveMesh.fullMesh,bc);
        end

        function createTractionSeparation(obj)
            s = obj.inputData;
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
            k.young   = ConstantFunction.create(obj.inputData.young,obj.cohesiveMesh.fullMesh);
            k.poisson = ConstantFunction.create(obj.inputData.poisson, obj.cohesiveMesh.fullMesh);
            k.cohesiveMesh = obj.cohesiveMesh;
            m  = Material.create(k);
            m.setCohesiveMesh(obj.cohesiveMesh);
        end

    end

end