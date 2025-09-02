classdef HardeningLawAT2 < HardeningLaw
 
    properties (Access = private)
        w
        qInf
        r1
    end

    methods (Access = public)
        
        function obj = HardeningLawAT2(cParams)
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
            qDeriv = obj.computeHardeningDerivative(r);
            qDot = (qDeriv.*(~isOverLimit));
        end
        
    end
    
    methods (Access = private)
        
        function initParams(obj,cParams)
            obj.w   = cParams.params.w;
            obj.r1   = cParams.params.r1;
            obj.qInf = obj.computeHardeningLimit();
        end

        function qInf = computeHardeningLimit(obj)
            qInf =  (2*obj.w*obj.r1 + obj.r1.^2-obj.r1^3)./(obj.w*2+obj.r1);
        end

        function itIs = isDamageOverLimit(obj,r)
            itIs = (r > obj.r1);
        end        

        function q = computeHardening(obj,r)
            q = (2*obj.w.*r)./(obj.w*2+r.^2);
        end

        function qDot = computeHardeningDerivative(obj,r)
            %qDot = (4*obj.w1^2 + 4*obj.w1*r - 8*obj.w1*r^2-r^4)./(2*obj.w1+r^2).^2;

            qDot = (2*obj.w.*(2*obj.w-r.^2))./(obj.w*2+r.^2).^2;
        end

    end
    
end