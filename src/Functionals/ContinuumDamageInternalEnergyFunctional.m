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
        end

        function sig = computeStress(obj,u,r)
            C   = obj.material.obtainTensorSecant(r);
            eps = SymGrad(u);
            sig = DDP(eps,C);
        end

        function [energy] = computeFunction(obj,u,r)
            sig = obj.computeStress(u,r);
            eps = SymGrad(u);
            en  = DDP(sig,eps);
            int = Integrator.compute(en,obj.mesh,obj.quadOrder);
            energy = 0.5*int;
        end

        function res = computeResidual(obj,u,r)
            stress = obj.computeStress(u,r);
            res = IntegrateRHS(@(v) DDP(SymGrad(v),stress),obj.test,obj.mesh,'Domain',obj.quadOrder);
        end

        function [Ktan,Ksec] = computeDerivativeResidual(obj,u,r)
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

        function LHS = computeLHS(obj,mat)
            LHS = IntegrateLHS(@(u,v) DDP(SymGrad(v),DDP(mat,SymGrad(u))),obj.test,obj.test,obj.mesh,'Domain',obj.quadOrder);
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