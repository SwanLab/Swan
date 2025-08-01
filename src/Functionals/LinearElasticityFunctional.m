classdef LinearElasticityFunctional < handle

    properties (Access = private)
        mesh
        material
    end

    methods (Access = public)

        function obj = LinearElasticityFunctional(cParams)
            obj.init(cParams)
        end

        function energy = computeCost(obj, uFun)
            C = obj.material;
            eps = SymGrad(uFun);
            fun = DDP(DDP(eps,C),eps);
            quadOrder = 3;
            energy = 0.5*(Integrator.compute(fun,obj.mesh,quadOrder));
        end

        function Ju = computeGradient(obj, uFun)
            eps = SymGrad(uFun);
            sig = DDP(obj.material,eps);
            test = LagrangianFunction.create(obj.mesh, uFun.ndimf, uFun.order);
            s.mesh = obj.mesh;
            s.quadratureOrder = 3;
            s.type = 'ShapeSymmetricDerivative';
            RHS = RHSIntegrator.create(s);
            Ju = RHS.compute(sig,test);
        end

        function Huu = computeHessian(obj, uFun)
            s.material = obj.material;
            s.test     = uFun;
            s.trial    = uFun;
            s.mesh     = obj.mesh;
            s.quadratureOrder = 3;
            s.type     = 'ElasticStiffnessMatrix';
            LHS = LHSIntegrator.create(s);
            Huu = LHS.compute();
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh     = cParams.mesh;
            obj.material = cParams.material;
        end

    end

end