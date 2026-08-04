classdef ContinuumDamageInternalEnergyFunctional < handle

    properties (Access = private)
        bMat, sMat, tMat
        mesh
        quadOrder
        test
    end

    properties (Access = private)
        RHS
    end

    methods (Access = public)

        function obj = ContinuumDamageInternalEnergyFunctional(cParams)
            obj.init(cParams);
        end

        function sig = computeStress(obj,u,r)
            C   = obj.sMat(r);
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
            C   = obj.bMat;
            eps = SymGrad(u);
            tau = sqrt(DDP(DDP(eps,C),eps));
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh      = cParams.mesh;
            obj.quadOrder = cParams.quadOrder;
            obj.bMat      = cParams.bMat;
            obj.sMat      = cParams.sMat;
            obj.tMat      = cParams.tMat;
            obj.test      = cParams.test;
        end

        function LHS = computeLHS(obj,mat)
            LHS = IntegrateLHS(@(u,v) DDP(SymGrad(v),DDP(mat,SymGrad(u))),obj.test,obj.test,obj.mesh,'Domain',obj.quadOrder);
        end

        function sec = computeDerivativeResidualSecant(obj,r)
            mat = obj.sMat(r);
            sec = obj.computeLHS(mat);
        end
        
        function tan = computeDerivativeResidualTangent(obj,u,r) 
            mat = obj.tMat(u,r);
            tan = obj.computeLHS(mat);
        end

    end
end