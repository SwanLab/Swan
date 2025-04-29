classdef TutorialXXPhaseField < handle

    properties (Access = public)
        initialGuess
        output
    end

    properties (Access = private)
        mesh
        boundaryConditions
        material
        matType
        dissipation
        functional
    end

    methods (Access = public)

        function obj = TutorialXXPhaseField()
            obj.defineCase();
            obj.createInitialGuess();
            obj.createMaterialPhaseField();
            obj.createDissipationInterpolation();
            obj.createPhaseFieldFunctional()
            obj.solvePhaseFieldProblem()
        end

    end

    methods (Access = private)

        function defineCase(obj)
            s.type.mesh = '1Elem';
            s.N         = 10; %Only for 'nElem' type.mesh
            s.type.bc   = 'displacementTraction';
            s.bcValues  = [0:0.001:0.1];
            [obj.mesh, obj.boundaryConditions] = BenchmarkManager.create(s);
        end

        function createInitialGuess(obj)
            u   = LagrangianFunction.create(obj.mesh,2,'P1');
            phi = LagrangianFunction.create(obj.mesh,1,'P1');
            %phi = obj.setInitialDamage(phi);
            obj.initialGuess.phi = phi;
            obj.initialGuess.u = u;
        end

        function phi = setInitialDamage(obj,phi)
            isInMiddle = obj.mesh.coord(:,1)>=0.5 & obj.mesh.coord(:,2)==0.5;
            fValues = phi.fValues;
            fValues(isInMiddle) = 0.01;
            phi.setFValues(fValues);
        end

        function createPhaseFieldFunctional(obj)
            s.mesh          = obj.mesh;
            s.material      = obj.material;
            s.dissipation   = obj.dissipation;
            s.l0            = 0.1;
            s.quadOrder     = 2;
            s.testSpace.u   = obj.initialGuess.u;
            s.testSpace.phi = obj.initialGuess.phi;
            s.energySplit   = isa(obj.material,'MaterialPhaseFieldAnalyticSplit');
            obj.functional  = PhaseFieldFunctional(s);
        end

        function createMaterialPhaseField(obj)
            E  = 210;
            nu = 0.3;

            s.type  = 'PhaseField';
            s.mesh  = obj.mesh;
            s.PFtype = 'Analytic';
            s.fileName = 'CircleMicroDamagePerimeter'; %Only for 'Homogenized' PFtype

            s.interp.interpolation = 'PhaseFieldDegradation';
            s.interp.degFunType    = 'AT';
            s.interp.ndim    = obj.mesh.ndim;
            s.interp.young   = ConstantFunction.create(E,obj.mesh);
            s.interp.poisson = ConstantFunction.create(nu,obj.mesh);

            obj.material = Material.create(s);
        end

        function createDissipationInterpolation(obj)
            s.mesh = obj.mesh;
            s.pExp = 2;
            Gc = 5e-3;
            if s.pExp == 1
                cw = 2/3;
            elseif s.pExp == 2
                cw = 1/2;
            end
            obj.dissipation.interpolation = PhaseFieldDissipationInterpolator(s);
            obj.dissipation.constant = Gc/(4*cw);
        end

        function outputData = solvePhaseFieldProblem(obj)
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
            outputData = PFComp.compute();
        end

    end

end