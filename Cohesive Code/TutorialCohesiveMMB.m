classdef TutorialCohesiveMMB < handle

    properties (Access = public)
        output  
        inputData
    end

    properties (Access = private)
        cohesiveFunctional
        solverType
        cohesiveMesh
        bc
    end

    methods (Access = public)

        function obj = TutorialCohesiveMMB()
            obj.init();
            obj.createMesh();
            [tractionSeparation] = obj.defineCase();
            obj.createCohesiveFunctional(tractionSeparation)
            obj.solveCohesiveProblem(tractionSeparation)
        end

        function solveCohesiveProblem(obj,ts)
            s.cohesiveMesh           = obj.cohesiveMesh;
            s.boundaryConditions     = obj.bc;
            s.functional             = obj.cohesiveFunctional;
            s.tolerance              = 1e-4;
            s.maxIter                = 100;
            s.tractionLaw            = ts;
            
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
            obj.inputData.young     = 120e3;
            obj.inputData.poisson   = 0.3;
    
            obj.inputData.lawType = 'TractionBiliniarCoupled';
            
            % TURON
            obj.inputData.tau0Normal       = 15;
            obj.inputData.tau0Shear        = 30;
            obj.inputData.firstCritEnergy  = 0.260;
            obj.inputData.secondCritEnergy = 1.002;
            obj.inputData.eta              = 2;
            obj.inputData.Kcoh        = 1e6;
            
            % obj.inputData.bcValues    = [0       :  5e-2    :  3.8,...
            %                              3.8     :  1e-3    :  10];

            obj.inputData.bcValues    = [0       :  5e-2    :  6,...
                                         6       :   2e-3   :  8,...
                                         8       :   2e-3   :  10];

            obj.inputData.problemType = 'MMB';

            obj.inputData.yCohLine = 1.55;
            obj.inputData.xCohLineMax = 150-35;

        end

        function createMesh(obj)
            a.fileName = 'B2GeometryRefined';
            t = FemDataContainer(a);
            s.baseMesh = t.mesh;

                s.xCohLineMax = obj.inputData.xCohLineMax; s.yCohLine = obj.inputData.yCohLine;
                obj.cohesiveMesh = NewCohesiveMesh(s);
        end

        function [ts] = defineCase(obj)
            obj.solverType = 'Newton';
            obj.createBoundaryConditions();
            ts = obj.createTractionSeparation();
        end

        function createBoundaryConditions(obj)
            boundaryConditions.type = obj.inputData.problemType;
            boundaryConditions.values = obj.inputData.bcValues;
            obj.bc  = BoundaryConditionsCreator(obj.cohesiveMesh.fullMesh,boundaryConditions);
        end

        function tractionSeparation = createTractionSeparation(obj)
            s = obj.inputData;
            tractionSeparation = CohesiveTractionSeparation(s);
        end

        function createCohesiveFunctional(obj,ts)
            s.tractionSeparation = ts;
            s.material           = obj.createMaterial();
            s.mesh               = obj.cohesiveMesh.fullMesh;
            s.cohesiveMesh       = obj.cohesiveMesh;
            s.quadOrder          = 2;
            s.test = LagrangianFunction.create(obj.cohesiveMesh.fullMesh,2,'P1');
            obj.cohesiveFunctional = CohesiveFunctional(s);
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