classdef ShFunc_ContinuumDamage < handle

   properties (Access = private)
        internalDamage
        internalElastic
        externalWork
        boundaryConditions
        quadOrder
    end

    methods (Access = public)

        function obj = ShFunc_ContinuumDamage(cParams)
            obj.createFunctionals(cParams);
        end

     %   function totalEnergy = computeEnergy(obj,quadOrder,u,r,fext)
     %       internalEnergy = obj.internalDamage.computeFunction(quadOrder,u,r);
     %       externalEnergy = obj.externalWork.computeFunction(u,fext,quadOrder);
     %       totalEnergy = internalEnergy - externalEnergy;
     %   end

        function res = computeResidual(obj,u,bc)
            fExt = bc.pointloadFun;
            Fext = obj.externalWork.computeGradient(u,fExt);
            Fint = obj.internalDamage.computeResidual(u,r);
            resT = Fint - Fext;
            res  = resT(bc.free_dofs);
        end

        function dRes = computeDerivativeResidual(obj,quadOrder,u,r)
            dRes = obj.internalDamage.computeDerivativeResidual(quadOrder,u,r);
        end

        function r = computeDamageEvolutionParam(obj,u)
            r = obj.internalDamage.computeDamageEvolutionParam(u);
        end

        function setROld(obj,r)
            obj.internalDamage.setROld(r);
        end

        function d = computeDamage(obj,r)
            d = obj.internalDamage.computeDamage(r);
        end
   end

    methods (Access = private)

        function createFunctionals(obj,cParams)
            obj.boundaryConditions = cParams.boundaryConditions;
            obj.internalDamage = ShFunc_ElasticDamage(cParams);
            %obj.internalElastic = ShFunc_Elastic(s);
            obj.externalWork = ShFunc_ExternalWork2(cParams);
        end

    end
end
