classdef ShFunc_ElasticDamage < handle
    
    properties (Access = private)
        material
        mesh
        r0
        H
    end
    
    properties (Access = private)
        rOld
    end
    
    methods (Access = public)
        
        function obj = ShFunc_ElasticDamage(cParams)
            obj.init(cParams)
        end
        
        function energy = computeFunction(obj,quadOrder,u,r)
            
            d = obj.computeDamage(r);

            Cdamage = obj.material.obtainTensor(d);
           
            epsi = SymGrad(u);
            funct = DDP(DDP(epsi,Cdamage),epsi);
            energy = 0.5*(Integrator.compute(funct,obj.mesh,quadOrder));
       
        end
        
        function res = computeResidual (obj,quadOrder,u,r)
            
            d = obj.computeDamage(r);
            Cdamage = obj.material.obtainTensor(d);

            epsi = SymGrad(u);
            b = DDP(epsi,Cdamage);
            
            S.type = 'ShapeSymmetricDerivative';
            S.quadratureOrder=quadOrder;
            S.mesh = obj.mesh;
            
            test = LagrangianFunction.create(obj.mesh, u.ndimf, u.order);
            
            rhs = RHSintegrator.create(S);
            
            res = rhs.compute(b,test);
            
        end
        
        function dRes = computeDerivativeResidual (obj,quadOrder,u,r)
            test = LagrangianFunction.create(obj.mesh, u.ndimf, u.order);
            d = obj.computeDamage(r);

            S.type = 'ElasticStiffnessMatrix';
            S.quadratureOrder = quadOrder;
            S.mesh = obj.mesh;
            S.material = obj.material.obtainTensor(d);
            S.test = test;
            S.trial = test;        
            
            lhs = LHSintegrator.create(S);            
            dRes = lhs.compute();
        end  
        
        function r = computeDamageEvolutionParam(obj,u)
            C = obj.material.obtainNonDamagedTensor;
            epsi = SymGrad(u);
            tauEpsi = power(DDP(DDP(epsi,C),epsi),0.5);
         
            if tauEpsi <= obj.rOld
                r = obj.rOld;
            else
                r = tauEpsi;
            end

        end

        function d = computeDamage(obj,r)
            q = obj.computeHardening();
            s.operation = @(xV) 1-(q.evaluate(r.evaluate(xV),obj.r0.evaluate(xV))./(r.evaluate(xV)));
            s.ndimf = 1;
            d = DomainFunction(s);

        end
        
        function setROld (obj,r)
            obj.rOld = r;        
        end
        
        function r = getROld (obj)
            r = obj.rOld;
        end
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.material = cParams.material;
            obj.mesh = cParams.mesh;
            obj.r0 = cParams.r0;
            obj.rOld = copy(obj.r0);
            obj.H = cParams.H;
        end

        function q = computeHardening(obj)
            s.operation = @(r,r0) r0 + obj.H *(r-r0);
            s.ndimf = 1;
            q = DomainFunction(s);
        end 
      
    end
    
end