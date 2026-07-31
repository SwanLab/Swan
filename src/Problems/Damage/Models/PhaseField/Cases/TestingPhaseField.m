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
        mat
        functional
    end

    methods (Access = public)

        function obj = TestingPhaseField(cParams)
            obj.init(cParams) 
            obj.defineCase();
            obj.createInitialGuess(cParams);
            obj.createMaterialPhaseField();
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
            fValues(:) = 0;
            phi.setFValues(fValues);
        end

        function phi = createDamageVariable(obj,phi)
            s.type = 'Damage';
            s.mesh = phi.mesh;
            s.fun  = phi;
            phi = DesignVariable.create(s);
        end

        function createPhaseFieldFunctional(obj)
            switch obj.matInfo.matType
                case 'Analytic'
                    s.energySplit = false;
                    s.C  = obj.mat.C;
                    s.dC = obj.mat.dC;
                    s.d2C = obj.mat.d2C;
                case 'AnalyticSplit'
                    s.energySplit = true;
                    s.mu  = obj.mat.mu;  s.dmu  = obj.mat.dmu;  s.d2mu  = obj.mat.d2mu;
                    s.k   = obj.mat.k;   s.dk   = obj.mat.dk;   s.d2k   = obj.mat.d2k;
                case 'Homogenized'
                    s.energySplit = false;
                    s.C  = obj.mat.C;
                    s.dC = obj.mat.dC;
                    s.d2C = obj.mat.d2C;
            end
            s.mesh          = obj.mesh;
            s.dissipation   = obj.createDissipationInterpolation();
            s.l0            = obj.l0;
            s.quadOrder     = 2;
            s.testSpace.u   = obj.initialGuess.u;
            s.testSpace.phi = obj.initialGuess.phi.fun;
            obj.functional  = PhaseFieldFunctional(s);
        end

        function createMaterialPhaseField(obj)
            N  = obj.mesh.ndim;
            E  = obj.matInfo.young;
            nu = obj.matInfo.poisson;

            switch obj.matInfo.matType
                case {'Analytic','AnalyticSplit'}
                    E0  = ConstantFunction.create(E,obj.mesh);
                    nu0 = ConstantFunction.create(nu,obj.mesh);
                    mu0 = LameParametersConverter.computeShearFromYoungAndPoisson(E0,nu0);
                    k0  = LameParametersConverter.computeBulkFromYoungAndPoisson(E0,nu0,N);

                    [mu,dmu,d2mu] = PhaseFieldInterpolators.compute(obj.matInfo.degradationType,mu0);
                    [k,dk,d2k]    = PhaseFieldInterpolators.compute(obj.matInfo.degradationType,k0);
                    l   = @(phi) LameParametersConverter.computeLambdaFromBulkAndShear(k(phi),mu(phi),N);
                    dl  = @(phi) LameParametersConverter.computeLambdaFromBulkAndShear(dk(phi),dmu(phi),N);
                    d2l = @(phi) LameParametersConverter.computeLambdaFromBulkAndShear(d2k(phi),d2mu(phi),N);

                    obj.mat.mu = mu; obj.mat.dmu = dmu; obj.mat.d2mu = d2mu;
                    obj.mat.k  = k;  obj.mat.dk  = dk;  obj.mat.d2k  = d2k;

                    muE   = @(phi) Expand(mu(phi),4);   lE   = @(phi) Expand(l(phi),4);
                    dmuE  = @(phi) Expand(dmu(phi),4);  dlE  = @(phi) Expand(dl(phi),4);
                    d2muE = @(phi) Expand(d2mu(phi),4); d2lE = @(phi) Expand(d2l(phi),4);

                    I   = ConstantFunction.create(eye4D(N),obj.mesh);
                    IxI = ConstantFunction.create(kronEye(N),obj.mesh);
                    obj.mat.C   = @(phi) 2*muE(phi).*I + lE(phi).*IxI;
                    obj.mat.dC  = @(phi) 2*dmuE(phi).*I + dlE(phi).*IxI;
                    obj.mat.d2C = @(phi) 2*d2muE(phi).*I + d2lE(phi).*IxI;

                case 'Homogenized'
                    s.fileName = obj.matInfo.fileName;
                    s.mesh     = obj.mesh;
                    s.young    = E;
                    hm = HomogenizedMaterialsReader(s);

                    obj.mat.C   = @(phi) hm.obtainTensor(phi);
                    obj.mat.dC  = @(phi) hm.obtainTensorDerivative(phi);
                    obj.mat.d2C = @(phi) hm.obtainTensorSecondDerivative(phi);
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