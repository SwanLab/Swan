classdef Tutorial07_1_PhaseFieldProblemAnalytical < handle

    properties (Access = public)
        degType
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
            switch obj.degType
                case 'AT'
                    s.energySplit = false;
                    s.C = obj.mat.C; s.dC = obj.mat.dC; s.d2C = obj.mat.d2C;
                case 'ATSplit'
                    s.energySplit = true;
                    s.mu = obj.mat.mu; s.dmu = obj.mat.dmu; s.d2mu = obj.mat.d2mu;
                    s.k = obj.mat.k; s.dk = obj.mat.dk; s.d2k = obj.mat.d2k;
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
            [mu,dmu,d2mu] = PhaseFieldInterpolators.compute(obj.degType,mu0);
            [k,dk,d2k]    = PhaseFieldInterpolators.compute(obj.degType,k0);
            l   = @(phi) k(phi) - (2/N)*mu(phi);
            dl  = @(phi) dk(phi) - (2/N)*dmu(phi);
            d2l = @(phi) d2k(phi) - (2/N)*d2mu(phi);
            obj.mat.mu = mu; obj.mat.dmu = dmu; obj.mat.d2mu = d2mu;
            obj.mat.k  = k; obj.mat.dk = dk; obj.mat.d2k = d2k;

            mu = @(phi) Expand(mu(phi),4); l = @(phi) Expand(l(phi),4);
            dmu = @(phi) Expand(dmu(phi),4); dl = @(phi) Expand(dl(phi),4);
            d2mu = @(phi) Expand(d2mu(phi),4); d2l = @(phi) Expand(d2l(phi),4);

            I           = ConstantFunction.create(eye4D(N),obj.mesh);
            IxI         = ConstantFunction.create(kronEye(N),obj.mesh);
            obj.mat.C   = @(phi) 2*mu(phi).*I + l(phi).*IxI;
            obj.mat.dC  = @(phi) 2*dmu(phi).*I + dl(phi).*IxI;
            obj.mat.d2C = @(phi) 2*d2mu(phi).*I + d2l(phi).*IxI;
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