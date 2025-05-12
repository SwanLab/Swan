classdef DamageLaw < handle

    properties (Access = private)
        hardeningLaw
    end
    
    methods (Access = public)

        function obj = DamageLaw(cParams)
            obj.init(cParams)
        end

        function d = computeFunction(obj,r,rOld)
            isDamaging = obj.checkDamaging(r,rOld);
            q = obj.hardeningLaw.computeFunction(r,isDamaging);
            d = (1-(q/r));
        end

        function dDot = computeDerivative(obj,r,rOld)
            isDamaging = obj.checkDamaging(r,rOld);
            q    = obj.hardeningLaw.computeFunction(r,isDamaging);
            qDot = obj.hardeningLaw.computeDerivative(r,isDamaging);
            dDot = (q - qDot*r)/(r.^3);
        end

        function qFun = getHardeningFun(obj,r)
            qFun = obj.hardeningLaw.computeFunction(r);
        end
        
    end
    
    methods (Access = private)

        function init(obj,cParams)
            obj.hardeningLaw = HardeningLaw.create(cParams);
        end

        function isDamaging = checkDamaging(~,r,rOld)
            isDamaging = (r > rOld);
        end
        
    end
end