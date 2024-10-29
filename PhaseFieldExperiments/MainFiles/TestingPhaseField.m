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
        functional
    end

    methods (Access = public)

        function obj = TestingPhaseField(cParams)
            obj.init(cParams) 
            obj.defineCase();
            obj.createInitialPhaseField();
            obj.createPhaseFieldFunctional()
        end

        function outputData = compute(obj)
            s.mesh = obj.mesh;
            s.boundaryConditions = obj.boundaryConditions;
            s.initPhi = obj.initialPhaseField;
            s.monitoring = obj.monitoring;
            s.functional = obj.functional;
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

        function createPhaseFieldFunctional(obj)
            s.mesh = obj.mesh;
            s.material = obj.createMaterialPhaseField();
            s.dissipation = obj.createDissipationInterpolation();
            s.l0 = obj.l0;
            s.quadOrder = 2;
            obj.functional = ShFunc_BrittlePhaseField(s);
        end

        function createInitialPhaseField(obj)
            phi = LagrangianFunction.create(obj.mesh,1,'P1');
            %phi.fValues(:) = 1e-12;
            obj.initialPhaseField = phi;
        end

        function material = createMaterialPhaseField(obj)
            s.mesh = obj.mesh;
            s.matInfo = obj.matInfo;
            s.Gc = obj.matInfo.Gc;
            s.type = obj.matInfo.matType;
            if s.type ~= "PhaseFieldAnalytic"
                s.fileName = obj.matInfo.fileName;
            end
            s.interpolation = obj.matInfo.degradation;
            material = Material.create(s);
        end

        function dissipation = createDissipationInterpolation(obj)
            s.pExp = obj.dissipInfo.pExp;
            dissipation.interpolation = PhaseFieldDissipationInterpolator(s);

            if s.pExp == 1
                dissipation.constant = obj.matInfo.Gc/(4*(2/3));
            elseif s.pExp == 2
                dissipation.constant = obj.matInfo.Gc/(4*(1/2));
            end
        end

    end

end