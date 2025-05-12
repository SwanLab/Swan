classdef ShFunc_ElasticDamage < handle

    properties (Access = private)
        material
        mesh
        quadOrder
    end

    properties (Access = private)
        RHS
        r
        rOld
        test
    end

    methods (Access = public)

        function obj = ShFunc_ElasticDamage(cParams)
            obj.init(cParams);
            obj.defineRHSIntegrator();
            obj.r    = LagrangianFunction.create(obj.mesh,1,'P0');
            obj.rOld = obj.defineROld(cParams.material.hardening.r0);
        end

        function setTestFunction(obj,u)
            obj.test = LagrangianFunction.create(obj.mesh, u.ndimf, u.order);
        end

        function [energy,C] = computeFunction(obj,u)
            C = obj.material.obtainTensorSecant(obj.r,obj.rOld);
            eps = SymGrad(u);
            sig = DDP(eps,C);
            en  = DDP(sig,eps);
            int = Integrator.compute(en,obj.mesh,obj.quadOrder);
            energy = 0.5*int;
        end

        function res = computeResidual(obj,u)
            C = obj.material.obtainTensorSecant(obj.r,obj.rOld);
            epsi = SymGrad(u);
            stress = DDP(epsi,C);
            res = obj.RHS.compute(stress,obj.test);
        end

        function [Ksec,Ktan] = computeDerivativeResidual(obj,u)
            Ksec = obj.computeDerivativeResidualSecant(obj.r,obj.rOld);
            Ktan = obj.computeDerivativeResidualTangent(u,obj.r,obj.rOld);
        end

        function updateDamageEvolutionParam(obj,u)
            C = obj.material.obtainNonDamagedTensor();
            epsi = SymGrad(u);
            tauEpsilon = sqrt(DDP(DDP(epsi,C),epsi));
            tauEpsilon = project(tauEpsilon,obj.r.order);

            fV = zeros(size(obj.r.fValues));

            nodesNoDamage = tauEpsilon.fValues <= obj.rOld.fValues;
            fV(nodesNoDamage) = obj.rOld.fValues(nodesNoDamage);
            fV(~nodesNoDamage) = tauEpsilon.fValues(~nodesNoDamage);

            obj.r.setFValues(fV);
        end

        %%%

        function setROld(obj)
            obj.rOld.setFValues(obj.r.fValues());
        end

        %%%

        function d = getDamage(obj)
            d = obj.material.getDamage();
        end

        function r = getR(obj)
            r = obj.r;
        end

        function qFun = getQ(obj)
            qFun = obj.material.getQ();
        end
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh      = cParams.mesh;
            obj.quadOrder = cParams.quadOrder;
            obj.material  = cParams.material;
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
        
        function defineRHSIntegrator(obj)
            s.type = 'ShapeSymmetricDerivative';
            s.quadratureOrder = obj.quadOrder;
            s.mesh  = obj.mesh;
            obj.RHS = RHSIntegrator.create(s);
        end

        function rOld = defineROld(obj,r0)
            rOld = copy(obj.r);
            rOld.setFValues(r0.*ones(size(rOld.fValues)));
        end

        function sec = computeDerivativeResidualSecant(obj)
            mat = obj.material.obtainTensorSecant(r,rOld);
            LHS = obj.createElasticLHS(mat);
            sec = LHS.compute();
        end
        
        function tan = computeDerivativeResidualTangent(obj,u,r,rOld) 
            isLoading = true;
            if (isLoading)
                C = obj.material.obtainTensorTanget(u,r,rOld); %Ctan
            else
                C = obj.material.obtainTensorSecant(r,rOld); %Csec
            end
            LHS = obj.createElasticLHS(C);
            tan = LHS.compute();
        end

    end
end