
classdef shFunc_ContinuumDamage < handle

    properties (GetAccess = public,SetAccess = private)
        r
    end
    properties (Access = private)
        material
        mesh

        internalDamage
        internalElastic
        externalWork
        
        H
    end

    methods (Access = public)

        function obj = shFunc_ContinuumDamage(cParams)
            obj.init(cParams);
            obj.createFunctionals();
            obj.setInitialDamage();
        end

        function totalEnergy = computeTotalEnergyDamage (obj, quadOrder, u,fext)
            internalEnergy = obj.internalDamage.computeFunction(quadOrder,u,obj.r);
            externalEnergy = obj.externalWork.computeFunction(u,fext,quadOrder);
            totalEnergy = internalEnergy - externalEnergy;
        end

        function totalEnergy = computeTotalEnergy (obj, quadOrder, u,fext)
            internalEnergy = obj.internalElastic.computeFunction(quadOrder,u);
            externalEnergy = obj.externalWork.computeFunction(u,fext,quadOrder);
            totalEnergy = internalEnergy - externalEnergy;
        end

        function F = computeJacobian(obj,quadOrder,u,bc,r)
            fExt = bc.pointloadFun;
            Fext = obj.externalWork.computeGradient(u,fExt,quadOrder);
            Fint = obj.internalDamage.computeJacobian(quadOrder,u,r);
            F = Fint - Fext;

        end

        function H = computeHessian(obj,quadOrder,u,r)
            H = obj.internalDamage.computeHessian(quadOrder,u,r);
        end

        function updateInternalVariableR(obj,uNew)
            obj.r = obj.internalDamage.updateDamage(obj.r,uNew);
        end
        function r = updateInternalVariableRAndReturnValue(obj,uNew)
            r = obj.internalDamage.updateDamage(obj.r,uNew);
        end
        
        function d = computeDamage (obj)
            d = obj.internalDamage.computeDamage(obj.r);
        end
        function setInternalVariableR(obj,r)
            obj.r = r;
        end


    end

    methods (Access = private)

        function init (obj,cParams)
            obj.material = cParams.material;
            obj.mesh = cParams.mesh;
            obj.r = cParams.r0;
            obj.H = cParams.H;
        end

        function createFunctionals(obj)
            s.mesh     = obj.mesh;
            s.material = obj.material;
            s.H = obj.H;
            s.r0 = obj.r;

            obj.internalDamage = shFunc_ElasticDamage(s);
            obj.internalElastic = shFunc_Elastic(s);
            obj.externalWork = shFunc_ExternalWork2(s);
        end
        function setInitialDamage(obj)
            obj.r = ConstantFunction.create(obj.r,obj.mesh);
        end
    end
end
