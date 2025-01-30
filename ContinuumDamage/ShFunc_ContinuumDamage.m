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

       % function totalEnergy = computeEnergy(obj,quadOrder,u,r,fext)
       %     internalEnergy = obj.internalDamage.computeFunction(quadOrder,u,r);
       %     externalEnergy = obj.externalWork.computeFunction(u,fext,quadOrder);
       %     totalEnergy = internalEnergy - externalEnergy;
       % end

        function [res] = computeResidual(obj,u,bc)
            fExt = bc.pointloadFun;
            Fext = obj.externalWork.computeResidual(u,fExt);
            Fint = obj.internalDamage.computeResidual(u);
            res = Fint - Fext;
            res  = res(bc.free_dofs);
        end

        function [K,dRes] = computeDerivativeResidual(obj,u,bc) 
            [K,dRes] = obj.internalDamage.computeDerivativeResidual(u);
            dRes = dRes(bc.free_dofs,bc.free_dofs);
        end

        function computeDamageEvolutionParam(obj,u)
            obj.internalDamage.computeDamageEvolutionParam(u);
        end

        function setROld(obj)
            obj.internalDamage.setROld();
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
