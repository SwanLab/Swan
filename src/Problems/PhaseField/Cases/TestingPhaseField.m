classdef TestingPhaseField < handle

    properties (Access = public)
        initialGuess
    end

    properties (Access = private)
        benchmark
        matInfo
        dissipInfo
        l0
        monitoring
        tolerance
        maxIter
        solver
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
            s.mesh               = obj.mesh;
            s.initialGuess       = obj.initialGuess;
            s.boundaryConditions = obj.boundaryConditions;
            s.functional         = obj.functional;
            s.monitoring         = obj.monitoring;
            s.tolerance          = obj.tolerance;
            s.maxIter            = obj.maxIter;
            s.solver             = obj.solver;
            PFComp = PhaseFieldComputer(s);

            outputData = PFComp.compute();
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.benchmark  = cParams.benchmark;
            obj.matInfo    = cParams.matInfo;
            obj.dissipInfo = cParams.dissipInfo;
            obj.l0         = cParams.l0;
            obj.monitoring = cParams.monitoring;
            obj.tolerance  = cParams.tolerance;
            obj.maxIter    = cParams.maxIter;
            obj.solver = cParams.solver;
        end

        function defineCase(obj)
            [obj.mesh, obj.boundaryConditions] = BenchmarkManager.create(obj.benchmark);
        end

        function createInitialGuess(obj,cParams)
            if isfield(cParams,'initialGuess')
                if isfield(cParams.initialGuess,'u')
                    u = cParams.initialGuess.u;
                else
                    u = LagrangianFunction.create(obj.mesh,2,'P1');
                end

                if isfield(cParams.initialGuess,'phi')
                    phi = cParams.initialGuess.phi;
                else
                    phi = LagrangianFunction.create(obj.mesh,1,'P1');
                    %phi = obj.setInitialDamage(phi);
                end
            else
                u = LagrangianFunction.create(obj.mesh,2,'P1');
                phi = LagrangianFunction.create(obj.mesh,1,'P1');
                %phi = obj.setInitialDamage(phi);
            end
            obj.initialGuess.u = u;
            obj.initialGuess.phi = obj.createDamageVariable(phi);
        end

        function phi = setInitialDamage(obj,phi)
            isInMiddle = obj.mesh.coord(:,1)>=0.5 & obj.mesh.coord(:,2)==0.5;
            fValues = phi.fValues;
            fValues(isInMiddle) = 0.01;
            %fValues = 0.01*ones(size(phi.fValues));
            phi.setFValues(fValues);
        end

        function phi = createDamageVariable(obj,phi)
            s.type = 'Damage';
            s.mesh = phi.mesh;
            s.fun  = phi;
            phi = DesignVariable.create(s);
        end

        function createPhaseFieldFunctional(obj)
            s.mesh          = obj.mesh;
            s.material      = obj.createMaterialPhaseField();
            s.dissipation   = obj.createDissipationInterpolation();
            s.l0            = obj.l0;
            s.quadOrder     = 2;
            s.testSpace.u   = obj.initialGuess.u;
            s.testSpace.phi = obj.initialGuess.phi.fun;
            s.energySplit   = (obj.matInfo.matType == "AnalyticSplit");
            obj.functional  = PhaseFieldFunctional(s);
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

    end

end