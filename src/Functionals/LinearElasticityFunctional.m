classdef LinearElasticityFunctional < handle

    properties (Access = private)
        mesh
        material
    end

    methods (Access = public)

        function obj = LinearElasticityFunctional(cParams)
            obj.init(cParams)
        end

        function energy = computeCost(obj,uFun,quadOrder)
            C = obj.material;
            eps = SymGrad(uFun);
            fun = DDP(DDP(eps,C),eps);
            energy = 0.5*(Integrator.compute(fun,obj.mesh,quadOrder));
        end

        function Ju = computeGradient(obj,uFun,quadOrder)
            C = obj.material;
            sigma = DDP(C,SymGrad(uFun));
            f = @(v) DDP(SymGrad(v),sigma);
            Ju = IntegrateRHS(f,uFun,obj.mesh,'Domain',quadOrder);
        end

        function Huu = computeHessian(obj,uFun,quadOrder)
            C     = obj.material;
            f = @(u,v) DDP(SymGrad(v),DDP(C,SymGrad(u)));
            Huu = IntegrateLHS(f,uFun,uFun,obj.mesh,'Domain',quadOrder);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh     = cParams.mesh;
            obj.material = cParams.material;
        end

    end

end