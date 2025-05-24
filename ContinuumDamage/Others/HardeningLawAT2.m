classdef HardeningLawAT2 < HardeningLaw
 
    properties (Access = private)
        w1
        qInf
        r1
    end

    methods (Access = public)
        
        function obj = HardeningLawAT2(cParams)
            obj@HardeningLaw(cParams)
            obj.initClassParams(cParams);
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
            qDeriv = obj.computeHardeningDerivative(r);
            qDot = (qDeriv.*(~isOverLimit));
        end
        
    end
    
    methods (Access = private)
        
        function initClassParams(obj,cParams)
            obj.w1        = cParams.w1;
            obj.r1       = cParams.r1;
            obj.qInf = obj.computeHardeningLimit();
        end

        function qInf = computeHardeningLimit(obj)
            qInf =  2*obj.w1*obj.r1/(obj.w1*2+obj.r1);
        end

        function itIs = isDamageOverLimit(obj,r)
            itIs = (r > obj.r1);
        end        

        function q = computeHardening(obj,r)
            q = 2*obj.w1*r/(obj.w1*2+r);
        end

        function qDot = computeHardeningDerivative(obj,r)
            qDot = (2*obj.w1*(2*obj.w1+r)- 2*obj.w1*r)/(2*obj.w1+r)^2;
        end

    end
    
end