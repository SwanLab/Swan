classdef HardeningLawLinear < HardeningLaw
 
    properties (Access = private)
        H
        r1
        qInf
    end

    methods (Access = public)
        
        function obj = HardeningLawLinear(cParams)
            obj@HardeningLaw(cParams)
            obj.initParams(cParams);
        end

        function qFun = computeFunction(obj,internalVariable)
            r = internalVariable.r;
            isOverLimit = obj.isDamageOverLimit(r);
            q = obj.computeHardening(r);
            qFun = (q.*(~isOverLimit) + obj.qInf.*isOverLimit);
        end

        function qDot = computeDerivative(obj,internalVariable)
            r = internalVariable.r;
            isOverLimit = obj.isDamageOverLimit(r);           
            qDot = (obj.H.*(~isOverLimit));
        end
        
    end
    
    methods (Access = private)
        
        function initParams(obj,cParams)
            obj.H    = cParams.params.hardening;
            obj.r1   = cParams.params.r1;
            obj.qInf = obj.computeHardeningLimit();
        end

        function qInf = computeHardeningLimit(obj)
            qInf = obj.r0 + obj.H*(obj.r1 - obj.r0);
        end

        function itIs = isDamageOverLimit(obj,r)
            itIs = (r > obj.r1);
        end        

        function q = computeHardening(obj,r)
            q = obj.r0 + obj.H*(r - obj.r0);
        end

    end
    
end