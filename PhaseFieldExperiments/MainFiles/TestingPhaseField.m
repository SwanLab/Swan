classdef TestingPhaseField < handle

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

        function obj = TestingPhaseField(cParams)
            obj.init(cParams) 
            obj.defineCase();
            obj.createInitialPhaseField();
            obj.createMaterialPhaseField();
            obj.createDissipationInterpolation();

        end

        function outputData = compute(obj)
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

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.monitoring = cParams.monitoring;
            obj.benchmark = cParams.benchmark;
            obj.matInfo = cParams.matInfo;
            obj.dissipInfo = cParams.dissipInfo;
            obj.l0 = cParams.l0;
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
            if s.type ~= "PhaseFieldAnalytic"
                s.fileName = obj.matInfo.fileName;
            end
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