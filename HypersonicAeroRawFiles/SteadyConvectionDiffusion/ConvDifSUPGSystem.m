classdef ConvDifSUPGSystem < handle

    properties (Access = private)
        mesh
        trial
        tau
    end

    methods (Access = public)
        function obj = ConvDifSUPGSystem(cParams)
            obj.init(cParams)
        end

        function [K,f] = compute(obj,a,nu,source)
            K = obj.computeLHS(a,nu);
            f = obj.computeRHS(a,source);
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh  = cParams.mesh;
            obj.trial = cParams.trial;
            obj.tau   = cParams.tau;
        end

        function K = computeLHS(obj,a,nu)
            Kadv = obj.computeAdvectionMatrix(a);
            Kst  = obj.computeStiffnessMatrix();
            Kast = obj.computeAnisotropicStiffnessMatrix();
            K    = Kadv + nu*Kst + a^2*Kast; % The a^2 can be seen as L2-norm squared
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

        function Kast = computeAnisotropicStiffnessMatrix(obj)
            s.trial = obj.trial;
            s.test  = obj.trial;
            s.mesh  = obj.mesh;
            s.type  = 'StabilizedStiffnessMatrix';
            s.quadratureOrder = 'QUADRATIC';
            s.Tau   = obj.tau;
            lhs     = LHSintegrator.create(s);
            Kast    = lhs.compute();
        end

        function f = computeRHS(obj,a,source)
            fShape = obj.computeRHSShape(source);
            fDeriv = obj.computeRHSShapeDeriv(source);
            t      = obj.tau;
            f      = fShape + t*a*fDeriv;
        end

        function f = computeRHSShape(obj,source)
            s.mesh     = obj.mesh;
            s.type     = 'ShapeFunction';
            s.quadType = 'QUADRATIC';
            int        = RHSintegrator.create(s);
            test       = obj.trial;
            f          = int.compute(source,test);
        end

        function f = computeRHSShapeDeriv(obj,source)
            s.mesh     = obj.mesh;
            s.type     = 'ShapeDerivative';
            s.quadType = 'QUADRATIC';
            int        = RHSintegrator.create(s);
            test       = obj.trial;
            f          = int.compute(source,test);
        end
    end
end