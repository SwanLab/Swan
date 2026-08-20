classdef Tutorial07_1_PhaseFieldProblemAnalytical < handle

    properties (Access = public)
        degType
        energySplit
        initialGuess
        output
    end

    properties (Access = private)
        mesh
        boundaryConditions
        mat
        dissipation
        functional
    end

    methods (Access = public)

        function obj = Tutorial07_1_PhaseFieldProblemAnalytical()
            obj.init()
            obj.defineCase();
            obj.createInitialGuess();
            obj.createMaterialPhaseField();
            obj.createDissipationInterpolation();
            obj.createPhaseFieldFunctional()
            obj.solvePhaseFieldProblem()
        end

    end

    methods (Access = private)

        function init(obj)
           close all;
           obj.degType = 'AT'; % 'ATSplit', ...
           obj.energySplit = false; %In the future may be different types, not just Vol-Dev
        end

        function defineCase(obj)
            s.mesh.type = '1Elem';
            s.bc.type   = 'DisplacementTractionY';
            s.bc.values  = [0:0.001:0.1];
            [obj.mesh, obj.boundaryConditions] = BenchmarkManager.create(s);
        end

        function createInitialGuess(obj)
            u   = LagrangianFunction.create(obj.mesh,2,'P1');
            phi = LagrangianFunction.create(obj.mesh,1,'P1');
            %phi = obj.setInitialDamage(phi);
            obj.initialGuess.phi = obj.createDamageVariable(phi);
            obj.initialGuess.u = u;
        end

        function phi = setInitialDamage(obj,phi)
            isInMiddle = obj.mesh.coord(:,1)>=0.5 & obj.mesh.coord(:,2)==0.5;
            fValues = phi.fValues;
            fValues(isInMiddle) = 0.01;
            phi.setFValues(fValues);
        end

        function phi = createDamageVariable(obj,phi)
            s.type = 'Damage';
            s.mesh = phi.mesh;
            s.fun  = phi;
            phi = DesignVariable.create(s);
        end

        function createPhaseFieldFunctional(obj)
            s.energySplit   = obj.energySplit;
            if obj.energySplit
                s.mu = obj.mat.mu; s.dmu = obj.mat.dmu; s.d2mu = obj.mat.d2mu;
                s.k = obj.mat.k; s.dk = obj.mat.dk; s.d2k = obj.mat.d2k;
            else
                s.C = obj.mat.C; s.dC = obj.mat.dC; s.d2C = obj.mat.d2C;
            end
            s.mesh          = obj.mesh;
            s.dissipation   = obj.dissipation;
            s.l0            = 0.1;
            s.quadOrder     = 2;
            s.testSpace.u   = obj.initialGuess.u;
            s.testSpace.phi = obj.initialGuess.phi.fun;
            obj.functional  = PhaseFieldFunctional(s);
        end

        function createMaterialPhaseField(obj)
            N   = obj.mesh.ndim;
            E0  = ConstantFunction.create(210,obj.mesh);
            nu0 = ConstantFunction.create(0.3,obj.mesh);
            mu0 = E0./(2*(1+nu0));
            k0  = E0./(N*(1-(N-1)*nu0));
            [mu,dmu,d2mu] = PhaseFieldDegradationInterpolator.compute(obj.degType,mu0);
            [k,dk,d2k]    = PhaseFieldDegradationInterpolator.compute(obj.degType,k0);
            obj.mat.mu = mu; obj.mat.dmu = dmu; obj.mat.d2mu = d2mu;
            obj.mat.k  = k; obj.mat.dk = dk; obj.mat.d2k = d2k;
            if obj.energySplit
                trc = @(u) squeezeParticular(trace(SymGrad(u)),1);
                trcSign = @(u) Heaviside(trc(u));

                obj.mat.k = @(phi,u) k(phi).*trcSign(u) + k0.*(1-trcSign(u));
                obj.mat.dk = @(phi,u) dk(phi).*trcSign(u);
                obj.mat.d2k = @(phi,u) d2k(phi).*trcSign(u);
            else
                l   = @(phi) LameParametersConverter.computeLambdaFromBulkAndShear(k(phi),mu(phi),N);
                dl  = @(phi) LameParametersConverter.computeLambdaFromBulkAndShear(dk(phi),dmu(phi),N);
                d2l = @(phi) LameParametersConverter.computeLambdaFromBulkAndShear(d2k(phi),d2mu(phi),N);

                muE   = @(phi) Expand(mu(phi),4);   lE   = @(phi) Expand(l(phi),4);
                dmuE  = @(phi) Expand(dmu(phi),4);  dlE  = @(phi) Expand(dl(phi),4);
                d2muE = @(phi) Expand(d2mu(phi),4); d2lE = @(phi) Expand(d2l(phi),4);

                I   = ConstantFunction.create(eye4D(N),obj.mesh);
                IxI = ConstantFunction.create(kronEye(N),obj.mesh);
                obj.mat.C   = @(phi) 2*muE(phi).*I + lE(phi).*IxI;
                obj.mat.dC  = @(phi) 2*dmuE(phi).*I + dlE(phi).*IxI;
                obj.mat.d2C = @(phi) 2*d2muE(phi).*I + d2lE(phi).*IxI;
            end
        end

        function createDissipationInterpolation(obj)
            s.mesh = obj.mesh;
            s.pExp = 2;
            Gc = 5e-3;
            if s.pExp == 1
                cw = 1/2;
            elseif s.pExp == 2
                cw = 2/3;
            end
            obj.dissipation.interpolation = PhaseFieldDissipationInterpolator(s);
            obj.dissipation.constant = Gc/(4*cw);
        end

        function solvePhaseFieldProblem(obj)
            s.mesh               = obj.mesh;
            s.initialGuess       = obj.initialGuess;
            s.boundaryConditions = obj.boundaryConditions;
            s.functional         = obj.functional;

            s.monitoring.set = true;
            s.monitoring.type = 'full';
            s.monitoring.print = true;

            s.tolerance.u = 1e-13;
            s.tolerance.phi = 1e-6;
            s.tolerance.stag = 1e-6;

            s.maxIter.u = 100;
            s.maxIter.phi = 300;
            s.maxIter.stag = 300;

            s.solver.type = 'Gradient';
            s.solver.tau  = 150;

            PFComp = PhaseFieldComputer(s);
            obj.output = PFComp.compute();
        end

    end

end