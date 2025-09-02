classdef HardeningLawAT1 < HardeningLaw
 
    properties (Access = private)
        w
        qInf
        r1
    end

    methods (Access = public)
        
        function obj = HardeningLawAT1(cParams)
            obj@HardeningLaw(cParams)
            obj.initParams(cParams);
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
        
        function initParams(obj,cParams)
            obj.w        = cParams.params.w;
            obj.r1       = cParams.params.r1;
            obj.qInf = obj.computeHardeningLimit();
        end

        function qInf = computeHardeningLimit(obj)
            qInf =  obj.w./(obj.r1);
        end

        function itIs = isDamageOverLimit(obj,r)
            itIs = (r > obj.r1);
        end     

        function itIs = isDamaging(obj,r)
            itIs = (r > obj.r0);
        end     

        function q = computeHardening(obj,r)
            q = obj.w./(r);
        end

        function qDot = computeHardeningDerivative(obj,r)
            qDot = -obj.w./r.^2;
        end

    end
    
end