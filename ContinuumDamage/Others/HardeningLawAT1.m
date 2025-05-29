classdef HardeningLawAT1 < HardeningLaw
 
    properties (Access = private)
        w1
        qInf
        r1
    end

    methods (Access = public)
        
        function obj = HardeningLawAT1(cParams)
            obj@HardeningLaw(cParams)
            obj.initClassParams(cParams);
        end

        function qFun = computeFunction(obj,internalVariable)
            r = internalVariable.r;
            isOverLimit = obj.isDamageOverLimit(r);
            isOverR0 = obj.isDamaging(r);
            q = obj.computeHardening(r);
            qFun = (~isOverR0).*obj.r0 + (isOverR0).*(q.*(~isOverLimit) + obj.qInf.*isOverLimit);
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
            qInf =  obj.w1./(obj.r1);
        end

        function itIs = isDamageOverLimit(obj,r)
            itIs = (r > obj.r1);
        end     

        function itIs = isDamaging(obj,r)
            itIs = (r > obj.r0);
        end     

        function q = computeHardening(obj,r)
            q = obj.w1./(r);
        end

        function qDot = computeHardeningDerivative(obj,r)
            qDot = -obj.w1./r.^2;
        end

    end
    
end