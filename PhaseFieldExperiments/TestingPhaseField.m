classdef TestingPhaseField < handle

    properties (Access = public)
        outputData
    end

    properties (Access = private)
        monitoring
        benchmark
        matInfo
        dissipInfo
        l0
    end

    properties (Access = private)
        mesh
        boundaryConditions
        initialPhaseField
        material
        dissipation
        constant
    end

    methods (Access = public)

        function obj = TestingPhaseField() %cParams
            obj.init() %cParams
            obj.defineCase();
            obj.createInitialPhaseField();
            obj.createMaterialPhaseField();
            obj.createDissipationInterpolation();
            obj.outputData = obj.solveProblem();
        end

    end

    methods (Access = private)

        function init(obj, ~)
            obj.monitoring.set = true;
            obj.monitoring.type = 'Reduced'; %'Full'
            obj.benchmark.type.mesh = '1Elem';
            obj.benchmark.N = 10;
            obj.benchmark.type.case = 'traction'; %'shear'
            obj.benchmark.type.bc = 'displacementTraction';
            obj.benchmark.bcValues = linspace(1e-4,1e-1,10);
            obj.matInfo.E  = 210;
            obj.matInfo.nu = 0.3;
            obj.matInfo.Gc = 5e-3;
            obj.matInfo.matType = 'PhaseFieldAnalytic'; %'PhaseFieldHomog'
            %obj.matInfo.fileName = 'IsoMicroDamage';
            obj.matInfo.degradation = 'PhaseFieldDegradation';
            obj.dissipInfo.type = 'PhaseFieldDissipationAT'; 
            obj.dissipInfo.pExp = 2;
            obj.l0 = 0.1;

            close all
            % obj.matInfo = cParams.matInfo;
            % obj.benchmark = cParams.benchmark;
            % obj.l0 = cParams.l0;
            % obj.pExp = cParams.pExp;
        end

        function outputData = solveProblem(obj)
            s.mesh = obj.mesh;
            s.initPhi = obj.initialPhaseField;
            s.material = obj.material;
            s.dissipation = obj.dissipation;
            s.boundaryConditions = obj.boundaryConditions;
            s.monitoring = obj.monitoring;
            s.l0 = obj.l0;
            PFComp = PhaseFieldComputer(s);

            outputData = PFComp.compute();
        end

        function defineCase(obj)
            [obj.mesh, obj.boundaryConditions] = benchmarkManager.create(obj.benchmark);
        end

        function createInitialPhaseField(obj)
            phi = LagrangianFunction.create(obj.mesh,1,'P1');
            obj.initialPhaseField = phi;
        end

        function createMaterialPhaseField(obj)
            s.mesh = obj.mesh;
            s.matInfo = obj.matInfo;
            s.Gc = obj.matInfo.Gc;
            s.type = obj.matInfo.matType;
            %s.fileName = obj.matInfo.fileName;
            s.interpolation = obj.matInfo.degradation;
            obj.material = Material.create(s);
        end

        function createDissipationInterpolation(obj)
            s.pExp = obj.dissipInfo.pExp;
            obj.dissipation.interpolation = PhaseFieldDissipationInterpolator(s);

            if s.pExp == 1
                obj.dissipation.constant = obj.matInfo.Gc/(4*(2/3));
            else%if s.pExp == 2
                obj.dissipation.constant = obj.matInfo.Gc/(4*(1/2));
            end
        end

    end

end