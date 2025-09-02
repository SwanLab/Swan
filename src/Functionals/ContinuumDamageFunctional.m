classdef ContinuumDamageFunctional < handle

    properties (Access = private)
        mesh
        quadOrder
        material
        test
    end

   properties (Access = private)
        internalEnergy
        externalWork
    end

    methods (Access = public)

        function obj = ContinuumDamageFunctional(cParams)
            obj.init(cParams);
            obj.createElasticDamage();
            obj.createExternalWork()
        end

        function totEnergy = computeEnergy(obj,u,r,bc)
           fExt = bc.pointloadFun;
           extE = obj.externalWork.computeCost(u,fExt,obj.quadOrder);
           intE = obj.internalEnergy.computeFunction(u,r);
           totEnergy = intE - extE;
       end

        function res = computeResidual(obj,u,r,bc)
            fExt = bc.pointloadFun;
            Fext = obj.externalWork.computeGradient(u,fExt,obj.quadOrder);
            Fint = obj.internalEnergy.computeResidual(u,r);
            res  = Fint - Fext;
        end

        function [Ktan,Ksec] = computeDerivativeResidual(obj,u,r) 
            [Ktan,Ksec] = obj.internalEnergy.computeDerivativeResidual(u,r);
        end

        function tau = computeTauEpsilon(obj,u)
            tau = obj.internalEnergy.computeTauEpsilon(u);
        end

        function q = getHardening(obj,r)
            q = obj.internalEnergy.getHardening(r);
        end

        function [d,sig] = getDamage(obj,u,r)
            d = obj.internalEnergy.getDamage(r);
            sig = obj.internalEnergy.computeStress(u,r);
        end
        
   end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh               = cParams.mesh;
            obj.quadOrder          = cParams.quadOrder;
            obj.material           = cParams.material;
            obj.test               = cParams.test;
        end

        function createElasticDamage(obj)
            s.mesh      = obj.mesh;
            s.quadOrder = obj.quadOrder;
            s.material  = obj.material;
            s.test      = obj.test;
            obj.internalEnergy = ContinuumDamageInternalEnergyFunctional(s);
        end

        function createExternalWork(obj)
            s.mesh        = obj.mesh;
            s.testSpace.u = obj.test;
            obj.externalWork = ExternalWorkFunctional(s);            
        end

    end
end
