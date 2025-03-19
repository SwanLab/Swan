classdef ShFunc_ContinuumDamage < handle

   properties (Access = private)
        internalDamage
        externalWork
    end

    methods (Access = public)

        function obj = ShFunc_ContinuumDamage(cParams)
            obj.createFunctionals(cParams);
        end

        function setTestFunctions(obj,u)           
            obj.internalDamage.setTestFunction(u);
            obj.externalWork.setTestFunction(u);
        end 

        function [totalEnergy,c] = computeEnergy(obj,u,bc)
           [internalEnergy,C] = obj.internalDamage.computeFunction(u);
           c = C.evaluate([0;0]);
           c = c(1,1);

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

        function [Ksec,dRes] = computeDerivativeResidual(obj,u,bc,control,index) 
            [K,Ksec] = obj.internalDamage.computeDerivativeResidual(u,control,index);
            dRes = K(bc.free_dofs,bc.free_dofs);
        end

        function computeDamageEvolutionParam(obj,u)
            obj.internalDamage.computeDamageEvolutionParam(u);
        end

        function setROld(obj)
            obj.internalDamage.setROld();
        end

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

        function createFunctionals(obj,cParams)
            obj.internalDamage = ShFunc_ElasticDamage(cParams);
            obj.externalWork = ShFunc_ExternalWork2(cParams);
        end

    end
end
