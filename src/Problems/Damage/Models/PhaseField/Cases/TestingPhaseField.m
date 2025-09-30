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
        initialDerivative
        functional
    end

    methods (Access = public)

        function obj = TestingPhaseField(cParams)
            obj.init(cParams) 
            obj.defineCase();
            obj.createInitialGuess(cParams);
            obj.computeInitialDerivative(cParams)
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
                    phi = obj.setInitialDamage(phi);
                end
            else
                u = LagrangianFunction.create(obj.mesh,2,'P1');
                phi = LagrangianFunction.create(obj.mesh,1,'P1');
                phi = obj.setInitialDamage(phi);
            end
            obj.initialGuess.u = u;
            obj.initialGuess.phi = obj.createDamageVariable(phi);
        end

        function phi = setInitialDamage(obj,phi)
            fValues = phi.fValues;
            fValues(:) = 1e-5;
            phi.setFValues(fValues);
        end

        function phi = createDamageVariable(obj,phi)
            s.type = 'Damage';
            s.mesh = phi.mesh;
            s.fun  = phi;
            phi = DesignVariable.create(s);
        end

        function computeInitialDerivative(obj,cParams)
            Gc = cParams.matInfo.Gc;
            E = cParams.matInfo.young;
            sigMax = 2;
            obj.initialDerivative = -2*(3/8)*(Gc/obj.l0)*E*(1/sigMax)^2;
        end

        function createPhaseFieldFunctional(obj)
            s.mesh          = obj.mesh;
            s.material      = obj.createMaterialPhaseField();
            s.dissipation   = obj.createDissipationInterpolation();
            s.l0            = obj.l0;
            s.quadOrder     = 3;
            s.testSpace.u   = obj.initialGuess.u;
            s.testSpace.phi = obj.initialGuess.phi.fun;
            s.energySplit   = (obj.matInfo.matType == "AnalyticSplit");
            obj.functional  = PhaseFieldFunctional(s);
        end

        function material = createMaterialPhaseField(obj)            
            s.type  = 'PhaseField';
            s.mesh  = obj.mesh;
            s.subType = obj.matInfo.matType;
            if s.subType == "Homogenized"
                s.fileName = obj.matInfo.fileName;
                s.young    = obj.matInfo.young;
            else
                s.interp = obj.defineDegradationFunction();  
            end
            material = Material.create(s);
        end

        function degParams = defineDegradationFunction(obj)
            E  = obj.matInfo.young;
            nu = obj.matInfo.poisson;
            ndim = obj.mesh.ndim;
            mu = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E,nu);
            kappa = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E,nu,ndim);

            degType = obj.matInfo.degradationType;
            degParams.interpolation = degType;
            switch degType
                case 'SIMPALL'
                    degParams.dim        = ndim;
                    degParams.matA.shear = 1e-5;
                    degParams.matA.bulk  = 1e-5;
                    degParams.matB.shear = mu;
                    degParams.matB.bulk  = kappa;
                case 'PhaseField'
                    if isfield(obj.matInfo,'degradationSubType')
                        degParams.subType    = obj.matInfo.degradationSubType;
                    end
                    degParams.ndim    = ndim;
                    degParams.young   = ConstantFunction.create(E,obj.mesh);
                    degParams.poisson = ConstantFunction.create(nu,obj.mesh);
                    degParams.initialDerivative = obj.initialDerivative;
            end
        end

        function dissipation = createDissipationInterpolation(obj)
            s.pExp = obj.dissipInfo.pExp;
            s.mesh = obj.mesh;
            dissipation.interpolation = PhaseFieldDissipationInterpolator(s);

            if s.pExp == 1
                dissipation.constant = (3/8)*obj.matInfo.Gc;
            elseif s.pExp == 2
                dissipation.constant = (1/2)*obj.matInfo.Gc;
            end
        end

    end

end