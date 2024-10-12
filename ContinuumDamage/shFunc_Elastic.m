classdef shFunc_Elastic < handle
    
    properties (Access = public)

    end
    
    properties (Access = private)
        material
        displacement
        mesh
    end
    
    methods (Access = public)
        
        function obj = shFunc_Elastic(cParams)
            obj.init(cParams)
            
        end
        
        function energy = computeFunction(obj,quadOrder)
            
            C = obj.material;
            epsi = SymGrad(obj.displacement);
            funct = DDP(DDP(epsi,C),epsi);
            energy = 0.5*(Integrator.compute(funct,obj.mesh,quadOrder));
       
        end
        
        function jacobian = computeJacobian(obj,quadOrder)
            
            C = obj.material;
            u = obj.displacement;
            
            epsi = SymGrad(u);
            b = DDP(epsi,C);
            
            
            S.type = 'ShapeSymmetricDerivative';
            S.quadratureOrder=1;
            S.mesh = obj.mesh;
            
            test = LagrangianFunction.create(obj.mesh, u.ndimf, u.order);
            
            rhs = RHSintegrator.create(S);
            
            jacobian = rhs.compute(b,test);
            
        end
        
        function hessian = computeHessian(obj,quadOrder)
            
            u = obj.displacement;
            test = LagrangianFunction.create(obj.mesh, u.ndimf, u.order);
            
            S.type = 'ElasticStiffnessMatrix';
            S.quadratureOrder = 2;
            S.mesh = obj.mesh;
            S.material = obj.material;
            S.test = test;
            S.trial = test;
           
            
            lhs = LHSintegrator.create(S);
            
            hessian = lhs.compute();
        
        end     
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.material = cParams.material;
            obj.displacement = cParams.u;
            obj.mesh = cParams.mesh;
        end
        
    end
    
end