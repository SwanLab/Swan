classdef shFunc_TotalEnergy < handle

    properties (Access = public)
    
    end
    properties (Access = private)
        material
        mesh

        internalDamage
        internalElastic
        external
        
        r0
        H
    end

    methods (Access = public)

        function obj = shFunc_TotalEnergy(cParams)
            obj.init(cParams);
            obj.createFunctionals ();
        end

        function totalEnergy = computeTotalEnergyDamage (obj, quadOrder, u,r,fext)
            internalEnergy = obj.internalDamage.computeFunction(quadOrder,u,r);
            externalEnergy = obj.external.computeFunction(u,fext,quadOrder);
            totalEnergy = internalEnergy + externalEnergy;
        end

        function totalEnergy = computeTotalEnergy (obj, quadOrder, u,fext)
            internalEnergy = obj.internalElastic.computeFunction(quadOrder,u);
            externalEnergy = obj.external.computeFunction(u,fext,quadOrder);
            totalEnergy = internalEnergy + externalEnergy;
        end

    end
    methods (Access = private)
        function init (obj,cParams)
            obj.material = cParams.material;
            obj.mesh = cParams.mesh;
            obj.r0 = cParams.r0;
            obj.H = cParams.H;
            
        end
        function createFunctionals(obj)
            s.mesh     = obj.mesh;
            s.material = obj.material;
            s.H = obj.H;
            s.r0 = obj.r0;

            obj.internalDamage = shFunc_ElasticDamage(s);
            obj.internalElastic = shFunc_Elastic(s)
            obj.external = shFunc_ExternalWork2(s);

        end
    end
end