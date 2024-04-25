classdef LHSintegratorConvDifSUSystem < handle

    properties (Access = private)
        mesh
        trial
    end

    methods (Access = public)
        function obj = LHSintegratorConvDifSUSystem(cParams)
            obj.init(cParams);
        end

        function K = compute(obj,a,nu)
            t    = obj.computeStabParameter(a,nu);
            Kadv = obj.computeAdvectionMatrix(a);
            Kst  = obj.computeStiffnessMatrix();
            Kast = obj.computeStabilizedStiffnessMatrix(t);
            K    = Kadv + nu*Kst + a^2*Kast;
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh  = cParams.mesh;
            obj.trial = cParams.trial;
        end

        function tau = computeStabParameter(obj,a,nu)
            h    = obj.mesh.computeMeanCellSize();
            Pe   = a*h/(2*nu);
            alfa = coth(Pe)-1/Pe;
            tau  = alfa*h/(2*a);
            if obj.trial.order == "P2"
                beta = 2*((coth(Pe)-1/Pe)-(cosh(Pe))^2*(coth(2*Pe)-1/(2*Pe)))/(2-(cosh(Pe))^2);
                tau_c = beta*h/(2*a);
                tau  = diag([tau_c,tau_c,tau]);
            end
        end

        function Kadv = computeAdvectionMatrix(obj,a)
            aFun = AnalyticalFunction.create(@(x) a*ones(size(x(1,:,:))),1,obj.mesh);
            s.trial    = obj.trial;
            s.test     = obj.trial;
            s.function = aFun;
            s.mesh     = obj.mesh;
            s.type     = 'AdvectionMatrixWithFunction';
            s.quadratureOrder = 'LINEAR';
            lhs        = LHSintegrator.create(s);
            Kadv       = lhs.compute();
        end

        function Kst = computeStiffnessMatrix(obj)
            s.trial = obj.trial;
            s.test  = obj.trial;
            s.mesh  = obj.mesh;
            s.type  = 'StiffnessMatrix';
            s.quadratureOrder = 'QUADRATIC';
            lhs     = LHSintegrator.create(s);
            Kst     = lhs.compute();
        end

        function Kast = computeStabilizedStiffnessMatrix(obj,tau)
            s.trial = obj.trial;
            s.test  = obj.trial;
            s.mesh  = obj.mesh;
            s.type  = 'StabilizedStiffnessMatrix';
            s.quadratureOrder = 'QUADRATIC';
            s.Tau   = tau;
            lhs     = LHSintegrator.create(s);
            Kast    = lhs.compute();
        end
    end
end