classdef ShFunc_ContinuumDamage < handle

    properties (Access = private)
        mesh
        quadOrder
        material
        test
    end

   properties (Access = private)
        internalDamage
        externalWork
    end

    methods (Access = public)

        function obj = ShFunc_ContinuumDamage(cParams)
            obj.init(cParams);
            obj.createElasticDamage();
            obj.createExternalWork()
        end

        function setTestFunctions(obj,u)           
            obj.internalDamage.setTestFunction(u);
            obj.externalWork.setTestFunction(u);
        end 

        function [totalEnergy] = computeEnergy(obj,u,r,bc)
           fExt = bc.pointloadFun;
           externalEnergy = obj.externalWork.computeFunction(u,fExt);
           [internalEnergy] = obj.internalDamage.computeFunction(u,r);
           totalEnergy = internalEnergy - externalEnergy;
       end

        function [res] = computeResidual(obj,u,r,bc)
            fExt = bc.pointloadFun;
            Fext = obj.externalWork.computeResidual(u,fExt);
            Fint = obj.internalDamage.computeResidual(u,r);
            res  = Fint - Fext;
            res  = res(bc.free_dofs);
        end

        function [dRes,Ksec] = computeDerivativeResidual(obj,u,r,bc) 
            [Ktan,Ksec] = obj.internalDamage.computeDerivativeResidual(u,r);
             % KsecFree = Ksec(bc.free_dofs,bc.free_dofs);
            dRes = Ktan(bc.free_dofs,bc.free_dofs);
        end

        function tau = computeTauEpsilon(obj,u)
            tau = obj.internalDamage.computeTauEpsilon(u);
        end

        function q = getHardening(obj,r)
            q = obj.internalDamage.getHardening(r);
        end

        function d = getDamage(obj,r)
            d = obj.internalDamage.getDamage(r);
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
            obj.internalDamage = ShFunc_ElasticDamage(s);
        end

        function createExternalWork(obj)
            s.mesh      = obj.mesh;
            s.quadOrder = obj.quadOrder;
            s.test      = obj.test;
            obj.externalWork = ShFunc_ExternalWork2(s);            
        end

    end
end
