classdef HardeningLawLinear < HardeningLaw
    properties (Access = private)
        H
        r1
    end

    methods (Access = public)
        
        function obj = HardeningLawLinear(cParams)
            obj@HardeningLaw(cParams)
            obj.initClassParams(cParams);
        end

        function q = computeQ(obj)
            qinf = qinfComputer();
            op = @(xV) obj.r0(xV) + obj.H*(obj.r(xV) - obj.r0(xV));

            q = @(xV) op(xV).*obj.isDamaging(xV)*~obj.isOverR1(xV) + qinf()*obj.isDamaging(xV)*obj.isOverR1(xV);
        end

        function qDot = computeQDerivative(obj)
            qDot = @(xV) obj.H.*obj.isDamaging(xV).*~obj.isOverR1(xV);
        end
        
    end
    
    methods (Access = private)
        
        function initClassParams(obj,cParams)
            obj.H = cParams.H;
            % obj.r1 = cParams.r1;
            obj.r1 = @(xV) cParams.r1.evaluate(xV);
        end

        function qinf = qinfComputer(obj)
            % qinf = @(xV) obj.r0.evaluate(xV) + obj.H*(r.evaluate(xV) - obj.r0.evaluate(xV));
            qinf = @(xV) obj.r0(xV) - obj.H*(obj.r1(xV) - obj.r0(xV));
        end
        
    end
    
end