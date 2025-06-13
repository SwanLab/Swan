classdef LinearElasticityFunctional < handle
    
    properties (Access = private)
        mesh
        material
    end
    
    methods (Access = public)
        
        function obj = LinearElasticityFunctional (cParams)
            obj.init(cParams)
        end

        function energy = compute(obj, uFun)
            C = obj.material;
            epsi = SymGrad(uFun);
            fun = DDP(DDP(epsi,C),epsi);
            quadOrder = 3;
            energy = 0.5*(Integrator.compute(fun,obj.mesh,quadOrder));
        end
        

        function Ju = computeGradient(obj, uFun)
            strain = SymGrad(uFun);
            sigma = DDP(obj.material,strain);
            test = LagrangianFunction.create(obj.mesh, uFun.ndimf, uFun.order);

            s.mesh = obj.mesh;
            s.quadratureOrder = 3;
            s.type = 'ShapeSymmetricDerivative';
            RHS = RHSIntegrator.create(s);
            Ju = RHS.compute(sigma,test);
        end

        function Huu = computeHessian(obj, uFun) 
            s.type     = 'ElasticStiffnessMatrix';
            s.mesh     = obj.mesh;
            s.material = obj.material;
            s.quadratureOrder = 3;
            s.test     = LagrangianFunction.create(obj.mesh,uFun.ndimf, 'P1');
            s.trial    = uFun;
            LHS = LHSIntegrator.create(s);
            Huu = LHS.compute();
        end

    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh   = cParams.mesh;
            obj.material = cParams.material;
        end

    end
    
end