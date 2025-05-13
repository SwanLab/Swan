classdef DamageLaw < handle

    properties (Access = private)
        internalDamageVariable
        hardeningLaw
    end
    
    methods (Access = public)

        function obj = DamageLaw(cParams)
            obj.init(cParams)
        end

        function d = computeFunction(obj,internalVariable)
            r = internalVariable.r;            
            q = obj.hardeningLaw.computeFunction(internalVariable);
            d = min(1-(q/r),(1-1e-10));
        end

        function dDot = computeDerivative(obj,internalVariable)
            r = internalVariable.r;
            isDamaging  = internalVariable.isDamaging();
            q    = obj.hardeningLaw.computeFunction(internalVariable);
            qDot = obj.hardeningLaw.computeDerivative(internalVariable);
            dDot = isDamaging.*(q - qDot.*r)./(r.^3);
        end

        function qFun = getHardening(obj,r)
            qFun = obj.hardeningLaw.computeFunction(r);
        end
        
    end
    
    methods (Access = private)

        function init(obj,cParams)
            obj.hardeningLaw = cParams.hardeningLaw;
        end
        
    end
end