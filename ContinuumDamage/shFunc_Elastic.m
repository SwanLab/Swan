classdef shFunc_Elastic < handle
    
    properties (Access = public)

    end
    
    properties (Access = private)
        material
        displacement
        mesh
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = shFunc_Elastic(cParams)
            obj.init(cParams)
            
        end
        
        function energy = computeFunction(obj,quadOrder)
            
            C = obj.material;
            epsi = symGrad(obj.displacement);
            funct = DDP(DDP(epsi,C),epsi);
            energy = 0.5*(Integrator.compute(funct,obj.mesh,quadOrder));
       
        end
        
        function jacobian = computeJacobian(obj,quadOrder)
            
            C = obj.material;
            epsi = symGrad(obj.displacement);
            b = DDP(epsi:C);
            u = obj.displacement;
            
            S.type = 'ShapeSymmetricDerivative';
            S.quadratureOrder=quadOrder;
            S.mesh = obj.mesh;
            test = LagrangianFunction.create(obj.mesh, u.ndimf, u.order);
            
            rhs = RHSintegrator.create(S);
            
            jacobian = rhs.compute(b,test);
            
        end
        
        function hessian = computeHessian(obj,quadOrder)
            
            test = LagrangianFunction.create(obj.mesh, u.ndimf, u.order);
            
            S.type = 'ElasticStiffnessMatrix';
            S.quadratureOrder = quadOrder;
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