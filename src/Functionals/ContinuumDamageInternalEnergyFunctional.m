classdef ContinuumDamageInternalEnergyFunctional < handle

    properties (Access = private)
        material
        mesh
        quadOrder
    end

    properties (Access = private)
        RHS
        test
    end

    methods (Access = public)

        function obj = ContinuumDamageInternalEnergyFunctional(cParams)
            obj.init(cParams);
            obj.createRHSIntegrator();
        end

        function [energy] = computeFunction(obj,u,r)
            C   = obj.material.obtainTensorSecant(r);
            eps = SymGrad(u);
            sig = DDP(eps,C);
            en  = DDP(sig,eps);
            int = Integrator.compute(en,obj.mesh,obj.quadOrder);
            energy = 0.5*int;
        end

        function sig = computeStress(obj,u,r)
            C   = obj.material.obtainTensorSecant(r);
            eps = SymGrad(u);
            sig = DDP(eps,C);
        end

        function res = computeResidual(obj,u,r)
            C      = obj.material.obtainTensorSecant(r);
            epsi   = SymGrad(u);
            stress = DDP(epsi,C);
            res = obj.RHS.compute(stress,obj.test);
        end

        function [Ktan,Ksec] = computeDerivativeResidual(obj,u,r)
            Ksec = obj.computeDerivativeResidualSecant(r);
            Ktan = obj.computeDerivativeResidualTangent(u,r);
        end

        function tau = computeTauEpsilon(obj,u)
            C = obj.material.obtainNonDamagedTensor();
            epsi = SymGrad(u);
            tau = sqrt(DDP(DDP(epsi,C),epsi));
            max(tau.evaluate([0;0]))
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

        function LHS = computeLHS(obj,mat)
            s.type = 'ElasticStiffnessMatrix';
            s.quadratureOrder = obj.quadOrder;
            s.mesh = obj.mesh;
            s.test  = obj.test;
            s.trial = obj.test;
            s.material = mat;
            integrator = LHSIntegrator.create(s);
            LHS = integrator.compute();
        end
        
        function createRHSIntegrator(obj)
            s.type = 'ShapeSymmetricDerivative';
            s.quadratureOrder = obj.quadOrder;
            s.mesh  = obj.mesh;
            obj.RHS = RHSIntegrator.create(s);
        end

        function sec = computeDerivativeResidualSecant(obj,r)
            mat = obj.material.obtainTensorSecant(r);
            sec = obj.computeLHS(mat);
        end
        
        function tan = computeDerivativeResidualTangent(obj,u,r) 
            mat = obj.material.obtainTensorTangent(u,r);
            tan = obj.computeLHS(mat);
        end

    end
end