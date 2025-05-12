classdef ShFunc_ContinuumDamage < handle

    properties (Access = private)
        mesh
        material
        quadOrder
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

        function [totalEnergy,c] = computeEnergy(obj,u,bc)
           [internalEnergy,C] = obj.internalDamage.computeFunction(u);
           %c = C.evaluate([0;0]);
           %c = c(1,1);

           fExt = bc.pointloadFun;
           externalEnergy = obj.externalWork.computeFunction(u,fExt);
           totalEnergy = internalEnergy - externalEnergy;
       end

        function [res] = computeResidual(obj,u,bc)
            fExt = bc.pointloadFun;
            Fext = obj.externalWork.computeResidual(u,fExt);
            Fint = obj.internalDamage.computeResidual(u);
            res = Fint - Fext;
            res  = res(bc.free_dofs);
        end

        function [dRes,Ksec] = computeDerivativeResidual(obj,u,bc,isLoading) 
            [K,Ksec] = obj.internalDamage.computeDerivativeResidual(u,isLoading);
            dRes = K(bc.free_dofs,bc.free_dofs);
        end

        function tau = computeTauEpsilon(obj,u)
            tau = obj.internalDamage.computeTauEpsilon(u);
        end

     %%%

        function setROld(obj)
            obj.internalDamage.setROld();
        end


     %%%


        function  r = getR(obj)
            r = obj.internalDamage.getR();
        end

        function q = getQ(obj)
            q = obj.internalDamage.getQ();
        end

        function d = getDamage(obj)
            d = obj.internalDamage.getDamage();
        end
   end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh               = cParams.mesh;
            obj.material           = cParams.material;
            obj.quadOrder          = cParams.quadOrder;
        end

        function createElasticDamage(obj)
            s.mesh      = obj.mesh;
            s.quadOrder = obj.quadOrder;
            s.material  = obj.material;            
            obj.internalDamage = ShFunc_ElasticDamage(s);
        end

        function createExternalWork(obj)
            s.mesh      = obj.mesh;
            s.quadOrder = obj.quadOrder;
            obj.externalWork = ShFunc_ExternalWork2(s);            
        end

    end
end
