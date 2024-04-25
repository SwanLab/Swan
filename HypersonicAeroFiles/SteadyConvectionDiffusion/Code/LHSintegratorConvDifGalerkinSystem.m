classdef LHSintegratorConvDifGalerkinSystem < handle

    properties (Access = private)
        mesh
        trial
    end

    methods (Access = public)
        function obj = LHSintegratorConvDifGalerkinSystem(cParams)
            obj.init(cParams)
        end

        function K = compute(obj,a,nu)
            Kadv = obj.computeAdvectionMatrix(a);
            Kst  = obj.computeStiffnessMatrix();
            K    = Kadv + nu*Kst;
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh  = cParams.mesh;
            obj.trial = cParams.trial;
        end

        function Kadv = computeAdvectionMatrix(obj,a)
            aFun       = AnalyticalFunction.create(@(x) a*ones(size(x(1,:,:))),1,obj.mesh);
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
    end
end