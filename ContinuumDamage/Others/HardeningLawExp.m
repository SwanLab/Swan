classdef HardeningLawExp < HardeningLaw     
    properties (Access = private)
        qInf
        A
    end 
    methods (Access = public)

        function obj = HardeningLawExp(cParams)
            obj@HardeningLaw(cParams)
            obj.initClassParams(cParams)  
        end

        function qFun = computeFunction(obj,internalVariable)
            r = internalVariable.r;
            qFun = obj.computeHardening(r);
        end

        function qDot = computeDerivative(obj,internalVariable)
            r = internalVariable.r;
            qDot = obj.computeHardeningDerivative(r);
        end
 
    end
    
    methods (Access = private)   
        function initClassParams(obj,cParams)
            obj.A = cParams.A;
            obj.qInf = cParams.qInf;
        end

        function q = computeHardening(obj,r)
            q = obj.qInf - (obj.qInf - obj.r0).*exp(obj.A*(1-(r./obj.r0)));
        end

        function qDot = computeHardeningDerivative(obj,r)
            qDot = obj.A*((obj.qInf - obj.r0)./(obj.r0)).*exp(obj.A*(1-(r/obj.r0)));
        end

    end
end