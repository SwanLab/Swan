classdef LinearElasticityFunctional < handle

    properties (Access = private)
        mesh
        material
        quadOrder
        test
    end

    methods (Access = public)

        function obj = LinearElasticityFunctional(cParams)
            obj.init(cParams)
        end

        function energy = computeCost(obj,uFun)
            C = obj.material;
            eps = SymGrad(uFun);
            fun = DDP(DDP(eps,C),eps);
            energy = 0.5*(Integrator.compute(fun,obj.mesh,obj.quadOrder));
        end

        function Ju = computeGradient(obj,uFun)
            C = obj.material;
            eps = SymGrad(uFun);
            Ju = IntegrateRHS(@(v) DDP(SymGrad(v),DDP(C,eps)), obj.test, obj.mesh,'Domain',obj.quadOrder);
        end

        function Huu = computeHessian(obj)
            C = obj.material;            
            Huu = IntegrateLHS(@(u,v) DDP(SymGrad(v),DDP(C,SymGrad(u))),obj.test,obj.test,obj.mesh,'Domain',obj.quadOrder);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh     = cParams.mesh;
            obj.material = cParams.material;
            obj.quadOrder = cParams.quadOrder;
            obj.test      = cParams.test;
        end

    end

end