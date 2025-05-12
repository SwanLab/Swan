classdef HardeningLawLinear < HardeningLaw
 
    properties (Access = private)
        H
        r1
        qInf
    end

    methods (Access = public)
        
        function obj = HardeningLawLinear(cParams)
            obj@HardeningLaw(cParams)
            obj.initClassParams(cParams);
        end

        function qFun = computeFunction(obj,internalVariable)
            r = internalVariable.r;
            isOverLimit = obj.isDamageOverLimit(r);
            isDamaging  = internalVariable.isDamaging();
            q = obj.computeHardening(internalVariable);
            qFun = isDamaging.*(q.*(~isOverLimit) + qinf.*isOverLimit);
        end

        function qDot = computeDerivative(obj)
            isOverLimit = obj.isDamageOverLimit();
            isDamaging  = internalVariable.isDamaging();            
            qDot = isDamaging.*(obj.H.*(~isOverLimit));
        end
        
    end
    
    methods (Access = private)
        
        function initClassParams(obj,cParams)
            obj.H        = cParams.H;
            obj.r1       = cParams.r1;
            obj.qInf = obj.computeHardeningLimit();
        end

        function qInf = computeHardeningLimit(obj)
            qInf = obj.r0 - obj.H*(obj.r1 - obj.r0);
        end

        function itIs = isDamageOverLimit(obj,r)
            itIs = (r > obj.r1);
        end        

        function q = computeHardening(obj,r)
            q = obj.r0 + obj.H*(r - obj.r0);
        end

    end
    
end