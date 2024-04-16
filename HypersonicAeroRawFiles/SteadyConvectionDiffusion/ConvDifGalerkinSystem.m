classdef ConvDifGalerkinSystem < handle

    properties (Access = private)
        mesh
        trial
    end

    methods (Access = public)
        function obj = ConvDifGalerkinSystem(cParams)
            obj.init(cParams)
        end

        function [K,f] = compute(obj,a,nu,source)
            K = obj.computeLHS(a,nu);
            f = obj.computeRHS(source);
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh  = cParams.mesh;
            obj.trial = cParams.trial;
        end

        function K = computeLHS(obj,a,nu)
            Kadv = obj.computeAdvectionMatrix(a);
            Kst  = obj.computeStiffnessMatrix();
            K    = Kadv + nu*Kst;
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

        function f = computeRHS(obj,source)
            s.mesh     = obj.mesh;
            s.type     = 'ShapeFunction';
            s.quadType = 'QUADRATIC';
            int        = RHSintegrator.create(s);
            test       = obj.trial;
            f          = int.compute(source,test);
        end
    end
end