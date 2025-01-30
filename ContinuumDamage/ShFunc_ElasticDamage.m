classdef ShFunc_ElasticDamage < handle
    
    properties (Access = private)
        material
        mesh
        test
        rhs
        r0
        H
        quadOrder
    end
    
    properties (Access = private)
        rOld
    end
    
    methods (Access = public)
        
        function obj = ShFunc_ElasticDamage(cParams)
            obj.init(cParams)
        end
        
        function energy = computeFunction(obj,u,r)            
            d = obj.computeDamage(r);
            C = obj.material.obtainTensor(d);           
            e  = SymGrad(u);
            s  = DDP(e,C);
            en = DDP(s,e);
            int = Integrator.compute(en,obj.mesh,obj.quadOrder);
            energy = 0.5*int;       
        end
        
        function res = computeResidual(obj,u,r)
            d = obj.computeDamage(r);
            C = obj.material.obtainTensor(d);
            e = SymGrad(u);
            s = DDP(e,C);
            res = obj.rhs.compute(s,obj.test);            
        end
        


        function dRes = computeDerivativeResidual (obj,quadOrder,u,r)
            S.d = obj.computeDamage(r);
            S.type = 'ElasticStiffnessMatrix';
            S.quadratureOrder = quadOrder;
            S.mesh = obj.mesh;
            S.test  = obj.test;
            S.trial = obj.test; 
            
            Sec = obj.computeSecandDResidual (S);
            Tan = obj.computeTangentDResidual (S,u,r);

            dRes = Sec% - Tan;
        end  
       
        function r = computeDamageEvolutionParam(obj,u) %PROBLEMA AQUI!
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
            s.operation = @(xV) 1-(q(r.evaluate(xV),obj.r0.evaluate(xV))./(r.evaluate(xV)));
            s.ndimf = 1;
            s.mesh  = obj.mesh;
            d = DomainFunction(s);
        end
        
        function setROld (obj,r)
            obj.rOld = r;        
        end
        
        function r = getROld(obj)
            r = obj.rOld;
        end

        function createTest(obj,u)           
            obj.test = LagrangianFunction.create(obj.mesh, u.ndimf, u.order);           
        end 
    end
   
    methods (Access = private)

        function init(obj,cParams)
            obj.quadOrder = 2;
            obj.material = cParams.material;
            obj.mesh = cParams.mesh;
            obj.r0 = cParams.r0;
            obj.rOld = copy(obj.r0);
            obj.H = cParams.H;
            obj.createRHSintegrator;
        end

        function createRHSintegrator(obj)
            s.type = 'ShapeSymmetricDerivative';
            s.quadratureOrder= obj.quadOrder;
            s.mesh = obj.mesh;             
            obj.rhs = RHSintegrator.create(s);
        end

        function q = computeHardening(obj)
            q = @(r,r0) r0 + obj.H *(r-r0);
        end 
        
        function Sec = computeSecandDResidual (obj,S)
            S.material = obj.material.obtainTensor(S.d);

            lhs = LHSintegrator.create(S);            
            Sec = lhs.compute();
        end 

        function Tan = computeTangentDResidual (obj,S,u,r)          
            q = obj.computeHardening;
            C = obj.material.obtainNonDamagedTensor;
            epsi = SymGrad(u);
            sig = DDP(epsi,C);
            
            S.material = C;
            S.fun = @(xV)((q(r.evaluate(xV),obj.r0.evaluate(xV))-obj.H*r.evaluate(xV))/((r.evaluate(xV))^3))*sig*sig;
            
            lhs = LHSintegrator.create(S);            
            Tan = lhs.compute();
        end 
    end    
end