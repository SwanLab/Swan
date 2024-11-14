classdef sh_TotalEnergy < handle

    properties (Access = public)
    
    end
    properties (Access = private)
        material
        mesh

        internal
        external
        
        r0
        H
    end

    methods (Access = public)
        function obj = sh_TotalEnergy(cParams)
            obj.init(cParams);
            obj.createFunctionals ();
        end
        function totalEnergy = computeTotalEnergy (obj, quadOrder, u,r,fext)
            internalEnergy = obj.internal.computeFunction(quadOrder,u,r);
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

            obj.internal = shFunc_ElasticDamage(s);
            obj.external = shFunc_ExternalWork2(s);

        end
    end
end