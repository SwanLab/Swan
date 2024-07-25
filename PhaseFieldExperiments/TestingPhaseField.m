classdef TestingPhaseField < handle

    properties (Access = private)
        benchmark;
        matInfo;
        pExp;
        l0;
    end

    properties (Access = private)
        mesh
        boundaryConditions
        initialPhaseField
        materialPhaseField
        dissipationPhaseField
        constant

        outputData
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
            obj.benchmark.type.mesh = '1Elem';
            obj.benchmark.type.case = 'traction'; %'shear'
            obj.benchmark.type.bc = 'displacementTraction';
            obj.benchmark.bcValues = linspace(1e-4,1e-1,100);
            obj.matInfo.E  = 210;
            obj.matInfo.nu = 0.3;
            obj.matInfo.Gc = 5e-3;
            obj.matInfo.matType = 'PhaseFieldAnalytic'; %'PhaseFieldHomog'
            obj.l0 = 0.1;
            obj.pExp = 2;

            close all
            % obj.matInfo = cParams.matInfo;
            % obj.benchmark = cParams.benchmark;
            % obj.l0 = cParams.l0;
            % obj.pExp = cParams.pExp;
        end

        function outputData = solveProblem(obj)
            s.mesh = obj.mesh;
            s.initialPhaseField = obj.initialPhaseField;
            s.materialPhaseField = obj.materialPhaseField;
            s.dissipationPhaseField = obj.dissipationPhaseField;
            s.boundaryConditions = obj.boundaryConditions;
            s.l0 = obj.l0;
            s.constant = obj.constant;
            outputData = PhaseFieldComputer(s);
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
            s.baseMaterial = obj.createBaseMaterial();
            s.degradation  = obj.createMaterialInterpolation();
            s.Gc = obj.matInfo.Gc;
            s.type = obj.matInfo.matType;
            obj.materialPhaseField = MaterialPhaseField(s);

            % sH.fileName    = 'IsoMicroDamage';
            % obj.materialPhaseField = HomogenizedPhaseField(sH);
            obj.materialPhaseField = Material.create(s);
        end

        function createDissipationInterpolation(obj)
            s.typeOfMaterial = 'ISOTROPIC';
            s.interpolation = 'PhaseFieldD';
            s.pExp = obj.pExp;
            obj.dissipationPhaseField = MaterialInterpolator.create(s);

            if s.pExp == 1
                obj.constant = obj.matInfo.Gc/(4*(2/3));
            else%if s.pExp == 2
                obj.constant = obj.matInfo.Gc/(4*(1/2));
            end
        end

        function matInt = createMaterialInterpolation(obj)
            s.typeOfMaterial = 'ISOTROPIC';
            s.interpolation = 'PhaseFieldI';
            matInt = MaterialInterpolator.create(s);
        end

    end

    methods (Access = private)

        function mat = createBaseMaterial(obj)
            sParam.mesh = obj.mesh;
            sParam.order = 'P1';
            sParam.fValues = ones(obj.mesh.nnodes,1)*obj.matInfo.E;
            young = LagrangianFunction(sParam);
            sParam.fValues = ones(obj.mesh.nnodes,1)*obj.matInfo.nu;
            poisson = LagrangianFunction(sParam);

            sIso.ndim = obj.mesh.ndim;
            sIso.young = young;
            sIso.poisson = poisson;
            mat = Isotropic2dElasticMaterial(sIso);
        end
    end

end