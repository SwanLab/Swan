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
            obj.createRHSIntegrator();
        end

        function setTestFunction(obj,u)           
            obj.test = LagrangianFunction.create(obj.mesh, u.ndimf, u.order);           
        end 
        
        function [energy,C] = computeFunction(obj,u)            
            obj.computeDamage();
            C = obj.material.obtainTensor(obj.d);           
            e  = SymGrad(u);
            s  = DDP(e,C);
            en = DDP(s,e);
            int = Integrator.compute(en,obj.mesh,obj.quadOrder);
            energy = 0.5*int;       
        end
        
        function res = computeResidual(obj,u)
            obj.computeDamage();
            C = obj.material.obtainTensor(obj.d);
            epsi = SymGrad(u);
            stress = DDP(epsi,C);
            res = obj.RHS.compute(stress,obj.test);            
        end
        
        function [K,Ksec] = computeDerivativeResidual(obj,u,control,index)
            obj.computeDamage();
            Ksec = obj.computeDerivativeResidualSecant();
            Ktan = obj.computeDerivativeResidualTangent(u,control,index);
            K = Ktan;
        end  
       
        function computeDamageEvolutionParam(obj,u)
            C = obj.material.obtainNonDamagedTensor();
            epsi = SymGrad(u);
            tauEpsilon = sqrt(DDP(DDP(epsi,C),epsi));
            tauEpsilon = project(tauEpsilon,obj.r.order);

            % tauEpsilon.evaluate([1;1])
            % epsi.evaluate([1;1])

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

        function createRHSIntegrator(obj)
            s.type = 'ShapeSymmetricDerivative';
            s.quadratureOrder= obj.quadOrder;
            s.mesh = obj.mesh;             
            obj.RHS = RHSIntegrator.create(s);
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
        
        function sec = computeDerivativeResidualSecant(obj)
            mat = obj.material.obtainTensor(obj.d);
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
            LHS = LHSIntegrator.create(s);
        end

        function tan = computeDerivativeResidualTangent(obj,u,control,index)    
            Csec = obj.material.obtainTensor(obj.d);
            if (index<control)
                C = obj.material.obtainNonDamagedTensor();
                epsi = SymGrad(u);
                sigBar = DDP(epsi,C);
                q = obj.computeHardening();
                op = @(xV)((q(obj.r.evaluate(xV),obj.r0.evaluate(xV))-obj.H*obj.r.evaluate(xV))./((obj.r.evaluate(xV)).^3));
                d_dot = DomainFunction.create(op,obj.mesh);
                Ctan2 = Expand(d_dot).*OP(sigBar,sigBar);
    
                op = @(xV) Csec.evaluate(xV) - Ctan2.evaluate(xV);
                Ctan = DomainFunction.create(op,obj.mesh);
            else
                Ctan = Csec;
            end
            LHS = obj.createElasticLHS(Ctan);         
            tan = LHS.compute();
        end 
    end    
end