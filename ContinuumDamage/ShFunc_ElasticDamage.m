classdef ShFunc_ElasticDamage < handle
    
    properties (Access = private)
        material
        mesh
        H
        r0
        quadOrder
    end
    
    properties (Access = private)
        test
        RHS

        rOld
        r
        d
    end
    
    methods (Access = public)

        function obj = ShFunc_ElasticDamage(cParams)
            obj.init(cParams)
            obj.rOld = copy(obj.r0);
            obj.r    = copy(obj.r0);
            obj.createRHSintegrator();
        end

        function setTestFunction(obj,u)           
            obj.test = LagrangianFunction.create(obj.mesh, u.ndimf, u.order);           
        end 
        
        % function energy = computeFunction(obj,u)            
        %     obj.computeDamage();
        %     C = obj.material.obtainTensor(obj.d);           
        %     e  = SymGrad(u);
        %     s  = DDP(e,C);
        %     en = DDP(s,e);
        %     int = Integrator.compute(en,obj.mesh,obj.quadOrder);
        %     energy = 0.5*int;       
        % end
        
        function res = computeResidual(obj,u)
            obj.computeDamage();
            C = obj.material.obtainTensor(obj.d);
            strain = SymGrad(u);
            stress = DDP(strain,C);
            res = obj.RHS.compute(stress,obj.test);            
        end
        
        function [K,dRes] = computeDerivativeResidual(obj,u)
            obj.computeDamage();
            Ksec = obj.computeDerivativeResidualSecant(obj.d);
            %Ktan = obj.computeDerivativeResidualTangent (s,u,obj.r);
            K = Ksec;
            dRes = Ksec;%Ktan
        end  
       
        function computeDamageEvolutionParam(obj,u)
            C = obj.material.obtainNonDamagedTensor;
            strain = SymGrad(u);
            tauEpsilon = power(DDP(DDP(strain,C),strain),0.5);
            tauEpsilon = project(tauEpsilon,obj.r.order);

            fV = zeros(size(obj.r.fValues));
            nodesNoDamage = tauEpsilon.fValues <= obj.rOld.fValues;
            fV(nodesNoDamage) = obj.rOld.fValues(nodesNoDamage);
            fV(~nodesNoDamage) = tauEpsilon.fValues(~nodesNoDamage);
            obj.r.setFValues(fV);
        end
        
        function setROld (obj)
            obj.rOld = obj.r;        
        end
        
        function d = getDamage(obj)
            d = obj.d;
        end
    end
   
    methods (Access = private)

        function init(obj,cParams)
            obj.material = cParams.material;
            obj.mesh = cParams.mesh;
            obj.quadOrder = cParams.quadOrder;
            obj.H = cParams.H;
            obj.r0 = cParams.r0;
        end

        function createRHSintegrator(obj)
            s.type = 'ShapeSymmetricDerivative';
            s.quadratureOrder= obj.quadOrder;
            s.mesh = obj.mesh;             
            obj.RHS = RHSintegrator.create(s);
        end

        function computeDamage(obj)
            q = obj.computeHardening();
            s.operation = @(xV) 1-(q(obj.r.evaluate(xV),obj.r0.evaluate(xV))./(obj.r.evaluate(xV)));
            s.ndimf = 1;
            s.mesh  = obj.mesh;
            obj.d = DomainFunction(s);
        end

        function q = computeHardening(obj)
            q = @(r,r0) r0 + obj.H *(r-r0);
        end 
        
        function sec = computeDerivativeResidualSecant(obj,d)
            mat = obj.material.obtainTensor(d);
            LHS = obj.createElasticLHS(mat);         
            sec = LHS.compute();
        end 

        function LHS = createElasticLHS(obj,material)
            s.type = 'ElasticStiffnessMatrix';
            s.quadratureOrder = obj.quadOrder;
            s.mesh = obj.mesh;
            s.test  = obj.test;
            s.trial = obj.test;
            s.material = material;
            LHS = LHSintegrator.create(s);
        end

        % function Tan = computeDerivativeResidualTangent (obj,s,u,r)          
        %     q = obj.computeHardening;
        %     C = obj.material.obtainNonDamagedTensor;
        %     epsi = SymGrad(u);
        %     sig = DDP(epsi,C);
        % 
        %     s.material = C;
        %     s.fun = @(xV)((q(r.evaluate(xV),obj.r0.evaluate(xV))-obj.H*r.evaluate(xV))/((r.evaluate(xV))^3))*sig*sig;
        % 
        %     lhs = LHSintegrator.create(s);            
        %     Tan = lhs.compute();
        % end 
    end    
end