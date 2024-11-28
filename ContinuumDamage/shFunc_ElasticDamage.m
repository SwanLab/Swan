classdef shFunc_ElasticDamage < handle
    
    properties (Access = public)

    end
    
    properties (Access = private)
        material
        mesh
        
        r0
        H

    end
    
    methods (Access = public)
        
        function obj = shFunc_ElasticDamage(cParams)
            obj.init(cParams)
            
        end
        
        function energy = computeFunction(obj,quadOrder,u,r)
            
            d = obj.computeDamage(r);

            Cdamage = obj.material.obtainTensor(d);
           
            epsi = SymGrad(u);
            funct = DDP(DDP(epsi,Cdamage),epsi);
            energy = 0.5*(Integrator.compute(funct,obj.mesh,quadOrder));
       
        end
        
        function jacobian = computeJacobian(obj,quadOrder,u,r)
            
            d = obj.computeDamage(r);
            Cdamage = obj.material.obtainTensor(d);

            epsi = SymGrad(u);
            b = DDP(epsi,Cdamage);
            
            S.type = 'ShapeSymmetricDerivative';
            S.quadratureOrder=quadOrder;
            S.mesh = obj.mesh;
            
            test = LagrangianFunction.create(obj.mesh, u.ndimf, u.order);
            
            rhs = RHSintegrator.create(S);
            
            jacobian = rhs.compute(b,test);
            
        end
        
        function hessian = computeHessian(obj,quadOrder,u,r)
            
            test = LagrangianFunction.create(obj.mesh, u.ndimf, u.order);
            d = obj.computeDamage(r);

            S.type = 'ElasticStiffnessMatrix';
            S.quadratureOrder = 2;
            S.mesh = obj.mesh;
            S.material = obj.material.obtainTensor(d);
            S.test = test;
            S.trial = test;
           
            
            lhs = LHSintegrator.create(S);
            
            hessian = lhs.compute();
        
        end  
        
        function rOut = updateDamage (obj,rIn,u)
    
            C = obj.material.obtainNonDamagedTensor;
            epsi = SymGrad(u);
            tauEpsi = power(DDP(DDP(epsi,C),epsi),0.5);
         
            if tauEpsi <= rIn

                rOut = rIn;
            else
                rOut = tauEpsi;
            end

        end

        function d = computeDamage(obj,r)
            q = obj.computeHardening();
            s.operation = @(xV) 1-(q.evaluate(r.evaluate(xV))./(r.evaluate(xV)));
            s.ndimf = 1;
            d = DomainFunction(s);

        end
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.material = cParams.material;
            obj.mesh = cParams.mesh;
            obj.r0 = cParams.r0;
            obj.H = cParams.H;
        end

        function q = computeHardening(obj)
            
            s.operation = @(r) obj.r0 + obj.H *(r-obj.r0);
            s.ndimf = 1;
            q = DomainFunction(s);

        end 
      
    end
    
end