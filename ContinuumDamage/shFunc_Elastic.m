classdef shFunc_Elastic < handle
    
    properties (Access = private)
        material
        mesh
    end
    
    methods (Access = public)
        
        function obj = shFunc_Elastic(cParams)
            obj.init(cParams) 
        end
        
        function energy = computeFunction(obj,quadOrder,u)
            C = obj.material;
            epsi = SymGrad(u);
            fun = DDP(DDP(epsi,C),epsi);
            energy = 0.5*(Integrator.compute(fun,obj.mesh,quadOrder));
       
        end
        
        function jacobian = computeJacobian(obj,quadOrder,u)
            s.type = 'ShapeSymmetricDerivative';
            s.quadratureOrder = quadOrder;
            s.mesh = obj.mesh;
            rhs = RHSintegrator.create(s);

            C = obj.material;
            sigma = DDP(C,SymGrad(u));
            test = LagrangianFunction.create(obj.mesh, u.ndimf, u.order);
            jacobian = rhs.compute(sigma,test);
        end
        
        function hessian = computeHessian(obj,quadOrder,u)
            test = LagrangianFunction.create(obj.mesh, u.ndimf, u.order);
            
            s.type = 'ElasticStiffnessMatrix';
            s.quadratureOrder = quadOrder;
            s.mesh = obj.mesh;
            s.material = obj.material;
            s.test = test;
            s.trial = test;
            lhs = LHSintegrator.create(s);
            hessian = lhs.compute();
        end     
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.material = cParams.material;
            obj.mesh = cParams.mesh;
        end
        
    end
    
end