classdef ShFunc_ElasticDamage < handle

    properties (Access = private)
        material
        mesh
        quadOrder
    end

    properties (Access = private)
        LHS
        RHS
        test
    end

    methods (Access = public)

        function obj = ShFunc_ElasticDamage(cParams)
            obj.init(cParams);
            obj.createRHSIntegrator();
            obj.createLHSIntegrator();
        end

        function [energy] = computeFunction(obj,u,r)
            C   = obj.material.obtainTensorSecant(r);
            eps = SymGrad(u);
            sig = DDP(eps,C);
            en  = DDP(sig,eps);
            int = Integrator.compute(en,obj.mesh,obj.quadOrder);
            energy = 0.5*int;
        end

        function res = computeResidual(obj,u,r)
            C      = obj.material.obtainTensorSecant(r);
            epsi   = SymGrad(u);
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

        function qFun = getHardening(obj,r)
            qFun = obj.material.getHardening(r);
        end

        function d = getDamage(obj,r)
            d = obj.material.getDamage(r);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh      = cParams.mesh;
            obj.quadOrder = cParams.quadOrder;
            obj.material  = cParams.material;
            obj.test      = cParams.test;
        end

        function createLHSIntegrator(obj)
            s.type = 'ElasticStiffnessMatrix';
            s.quadratureOrder = obj.quadOrder;
            s.mesh = obj.mesh;
            s.test  = obj.test;
            s.trial = obj.test;
            obj.LHS = LHSIntegrator.create(s);
        end
        
        function createRHSIntegrator(obj)
            s.type = 'ShapeSymmetricDerivative';
            s.quadratureOrder = obj.quadOrder;
            s.mesh  = obj.mesh;
            obj.RHS = RHSIntegrator.create(s);
        end

        function sec = computeDerivativeResidualSecant(obj,r)
            mat = obj.material.obtainTensorSecant(r);
            sec = obj.LHS.compute(mat);
        end
        
        function tan = computeDerivativeResidualTangent(obj,u,r) 
            isLoading = true;
            if (isLoading)
                mat = obj.material.obtainTensorTangent(u,r); %Ctan
            else
                mat = obj.material.obtainTensorSecant(r); %Csec
            end
            tan = obj.LHS.compute(mat);
        end

    end
end