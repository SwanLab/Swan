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

        function res = computeResidual(obj,u,r)
            C = obj.material.obtainTensorSecant(r);
            epsi = SymGrad(u);
            stress = DDP(epsi,C);
            res = obj.RHS.compute(stress,obj.test);
        end

        function [Ksec,Ktan] = computeDerivativeResidual(obj,u,r)
            Ksec = obj.computeDerivativeResidualSecant(r);
            Ktan = obj.computeDerivativeResidualTangent(u,r);
        end

        function tau = computeTauEpsilon(obj,u)
            C = obj.material.obtainNonDamagedTensor();
            epsi = SymGrad(u);
            tau = sqrt(DDP(DDP(epsi,C),epsi));
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

        function sec = computeDerivativeResidualSecant(obj,r)
            mat = obj.material.obtainTensorSecant(r);
            LHS = obj.createElasticLHS(mat);
            sec = LHS.compute();
        end
        
        function tan = computeDerivativeResidualTangent(obj,u,r) 
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