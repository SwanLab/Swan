classdef ShFunc_ContinuumDamage < handle

   properties (Access = private)
        internalDamage
        internalElastic
        externalWork
        boundaryConditions
        quadOrder
        r
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
            Fext = obj.externalWork.computeGradient(u,fExt,obj.quadOrder);
            Fint = obj.internalDamage.computeResidual(u,obj.r);
            resT = Fint - Fext;
            res  = resT(bc.free_dofs); %THE FREE DOFS ARE RETURNED, revisar el pq
        end

        function dRes = computeDerivativeResidual(obj,u,bc) 
            dres = obj.internalDamage.computeDerivativeResidual(obj.quadOrder,u,obj.r);
            dRes = dres(bc.free_dofs,bc.free_dofs); %THE FREE DOFS ARE RETURNED, revisar el pq
        end

        function computeDamageEvolutionParam(obj,u)
            obj.r = obj.internalDamage.computeDamageEvolutionParam(u);
        end

        function setROld(obj)
            obj.internalDamage.setROld(obj.r);
        end

        function d = computeDamage(obj,r)
            d = obj.internalDamage.computeDamage(r);
        end

        function updateBoundaryConditions (obj,bc)
            obj.boundaryConditions = bc;
        end

        function createTest(obj,u)           
            obj.internalDamage.createTest(u);
        end 
   end

    methods (Access = private)

        function createFunctionals(obj,cParams)
            obj.boundaryConditions = cParams.boundaryConditions;
            obj.internalDamage = ShFunc_ElasticDamage(cParams);
            obj.quadOrder = cParams.quadOrder;
            %obj.internalElastic = ShFunc_Elastic(s);
            obj.externalWork = ShFunc_ExternalWork2(cParams);
        end

    end
end
