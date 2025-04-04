classdef TestingPhaseField < handle

    properties (Access = public)
        initialGuess
    end

    properties (Access = private)
        monitoring
        benchmark
        matInfo
        dissipInfo
        tolerance
        solverType
        l0
    end

    properties (Access = private)
        mesh
        boundaryConditions
        functional
    end

    methods (Access = public)

        function obj = TestingPhaseField(cParams)
            obj.init(cParams) 
            obj.defineCase();
            obj.createInitialGuess(cParams);
            obj.createPhaseFieldFunctional()
        end

        function outputData = compute(obj)
            s.mesh = obj.mesh;
            s.boundaryConditions = obj.boundaryConditions;
            s.initialGuess = obj.initialGuess;
            s.monitoring = obj.monitoring;
            s.functional = obj.functional;
            s.tolerance = obj.tolerance;
            s.solverType = obj.solverType;
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
            obj.tolerance = cParams.tolerance;
            obj.solverType = cParams.solverType;
            obj.l0 = cParams.l0;
        end

        function defineCase(obj)
            [obj.mesh, obj.boundaryConditions] = BenchmarkManager.create(obj.benchmark);
        end

        function createInitialGuess(obj,cParams)
            if isfield(cParams,'initialGuess')
                if isfield(cParams.initialGuess,'u')
                    obj.initialGuess.u = cParams.initialGuess.u;
                else
                    u = LagrangianFunction.create(obj.mesh,2,'P1');
                    obj.initialGuess.u = u;
                end

                if isfield(cParams.initialGuess,'phi')
                    obj.initialGuess.phi = cParams.initialGuess.phi;
                else
                    phi = LagrangianFunction.create(obj.mesh,1,'P1');
                    %phi = obj.setInitialDamage(phi);
                    obj.initialGuess.phi = phi;
                end
            else
                u = LagrangianFunction.create(obj.mesh,2,'P1');
                phi = LagrangianFunction.create(obj.mesh,1,'P1');
                %phi = obj.setInitialDamage(phi);
                obj.initialGuess.phi = phi;
                obj.initialGuess.u = u;
            end
        end

        function createPhaseFieldFunctional(obj)
            s.mesh          = obj.mesh;
            s.material      = obj.createMaterialPhaseField();
            s.dissipation   = obj.createDissipationInterpolation();
            s.l0            = obj.l0;
            s.quadOrder     = 2;
            s.testSpace.u   = obj.initialGuess.u;
            s.testSpace.phi = obj.initialGuess.phi;
            s.energySplit   = (obj.matInfo.matType == "AnalyticSplit");
            obj.functional = PhaseFieldFunctional(s);
        end

        function material = createMaterialPhaseField(obj)
            E  = obj.matInfo.young;
            nu = obj.matInfo.poisson;
            
            s.type  = 'PhaseField';
            s.mesh  = obj.mesh;
            s.PFtype = obj.matInfo.matType;
            if s.PFtype == "Homogenized"
                s.fileName = obj.matInfo.fileName;
            else
                s.interp.interpolation = 'PhaseFieldDegradation';
                s.interp.degFunType    = obj.matInfo.degradationType;
                s.interp.ndim    = obj.mesh.ndim;
                s.interp.young   = ConstantFunction.create(E,obj.mesh);
                s.interp.poisson = ConstantFunction.create(nu,obj.mesh);
            end
            material = Material.create(s);
        end

        function dissipation = createDissipationInterpolation(obj)
            s.pExp = obj.dissipInfo.pExp;
            s.mesh = obj.mesh;
            dissipation.interpolation = PhaseFieldDissipationInterpolator(s);

            if s.pExp == 1
                dissipation.constant = obj.matInfo.Gc/(4*(2/3));
            elseif s.pExp == 2
                dissipation.constant = obj.matInfo.Gc/(4*(1/2));
            end
        end

        function phi = setInitialDamage(obj,phi)
            isInMiddle = obj.mesh.coord(:,1)>=0.5 & obj.mesh.coord(:,2)==0.5;
            fValues = phi.fValues;
            fValues(isInMiddle) = 0.01;
            %fValues = 0.01*ones(size(phi.fValues));
            phi.setFValues(fValues);
        end

    end

end