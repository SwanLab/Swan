classdef TutorialCohesive < handle

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

        function obj = TutorialCohesive()
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
            obj.inputData.young     = 120e20;
            obj.inputData.poisson   = 0;
    
            obj.inputData.lawType = 'TractionBiliniarCoupled';
            
            % TURON
            obj.inputData.tau0Normal       = 15e6;
            obj.inputData.tau0Shear        = 30e6;
            obj.inputData.firstCritEnergy  = 5*260;
            obj.inputData.secondCritEnergy = 1002;
            obj.inputData.eta              = 2;
                        
            % NO TURON
            obj.inputData.jumpCrit         = 1.25e-7;
            obj.inputData.jumpFinal        = 0.025e-3;

            % % UnitElem
            obj.inputData.Kcoh        = 1e12;

            % ======= UNLOAD ==========

            % 15
            % obj.inputData.bcValues    = [0       :  1e-6    :  5e-5,...
            %                             5e-5  :    -2e-6   :  0 ,...
            %                             0        :   -2e-6  :  -3e-5];

            % 30
            % obj.inputData.bcValues    = [0       :  1e-6    :  6e-5,...
            %                             6e-5  :    -2e-6   :  0 ,...
            %                             0        :   -2e-6  :  -1.7e-5];


            % ======= TOTAL ============

            % 15
            obj.inputData.bcValues    = [0       :  1e-6    :  6.7e-5];

            % 30
            % obj.inputData.bcValues    = [0       :  1e-6    :  8.6e-5];




            obj.inputData.problemType = 'SingleShear';
            obj.inputData.l = 1;
            obj.inputData.h = 1;
            obj.inputData.yCohLine = 0;
            obj.inputData.xCohLineMax =1000;
            obj.inputData.nx = 1;
            obj.inputData.ny = 1;

            % % DCB
            % obj.inputData.Kcoh      = 1e13;
            % obj.inputData.bcValues = [0:0.0001:1.2e-3,1.2e-3:0.00001:1.5e-3, 1.5e-3:0.00005:7e-3, 7e-3:0.00001:10e-3]; % Amb Gc = Gc
            % % obj.inputData.bcValues = [0:0.0001:2.2e-3,2.2e-3:0.00002:3.5e-3, 3.5e-3:0.00002:10e-3]; % Amb Gc = 10Gc
            % obj.inputData.problemType = 'DoubleCantileverBeam';
            % obj.inputData.l = 150e-3;
            % obj.inputData.h = 1.55*2e-3;
            % obj.inputData.yCohLine = 0.5*obj.inputData.h;
            % obj.inputData.xCohLineMax = 115e-3;
            % obj.inputData.nx = 2000;
            % obj.inputData.ny = 10;


            % % ENF
            % obj.inputData.Kcoh      = 1e13;
            % obj.inputData.bcValues = [ ...
            %         0     : 1e-4     :6e-3,...
            %         6e-3  : 2e-5     :6.5e-3,...
            %         6.5e-3  : 1e-6     :9e-3];
            % 
            % obj.inputData.problemType = 'EndNotchedFlex';
            % obj.inputData.l = 150e-3;
            % obj.inputData.h = 1.55*2e-3;
            % obj.inputData.yCohLine = 0.5*obj.inputData.h;
            % obj.inputData.xCohLineMax = 115e-3;
            % obj.inputData.nx = 1000;
            % obj.inputData.ny = 10;

        end
    
        function createMesh(obj)
            s.xCohLineMax = obj.inputData.xCohLineMax; s.yCohLine = obj.inputData.yCohLine;
            s.baseMesh = QuadMesh(obj.inputData.l ,obj.inputData.h, obj.inputData.nx, obj.inputData.ny);
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