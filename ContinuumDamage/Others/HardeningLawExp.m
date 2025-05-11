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
        function q = computeQ(obj)
            op = @(xV) obj.qInf - (obj.qInf - obj.r0(xV))*exp(obj.A*(1-obj.r(xV)/obj.r0(xV)));
            q = @(xV) op(xV).*obj.isDamaging(xV);
        end
        function qDot = computeQDerivative(obj)
            op = @(xV) obj.A * (obj.qInf - obj.r0(xV))/(obj.r0(xV)) * exp(obj.A*(1-obj.r(xV)/obj.r0(xV)));
            qDot = @(xV) op(xV).*obj.isDamaging(xV);
        end  
    end
    
    methods (Access = private)   
        function initClassParams(obj,cParams)
            obj.A = cParams.A;
            obj.qInf = cParams.qInf;
        end
    end
end