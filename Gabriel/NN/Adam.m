classdef Adam < SGD

    properties (Access = private)
        beta1
        beta2
        epsilon
        firstMoment
        secondMoment
        nSteps
    end

    methods (Access = public)

        function obj = Adam(s)
            obj@SGD(s);
            obj.beta1   = Adam.readOption(s,'beta1',0.9);
            obj.beta2   = Adam.readOption(s,'beta2',0.999);
            obj.epsilon = Adam.readOption(s,'epsilon',1e-8);
            obj.nSteps  = 0;
        end

    end

    methods (Access = protected)

        function x = step(obj,x,alpha,grad)
            if isempty(obj.firstMoment)
                obj.firstMoment  = zeros(size(grad));
                obj.secondMoment = zeros(size(grad));
            end
            obj.nSteps       = obj.nSteps + 1;
            obj.firstMoment  = obj.beta1*obj.firstMoment  + (1-obj.beta1)*grad;
            obj.secondMoment = obj.beta2*obj.secondMoment + (1-obj.beta2)*grad.^2;
            mHat = obj.firstMoment /(1-obj.beta1^obj.nSteps);
            vHat = obj.secondMoment/(1-obj.beta2^obj.nSteps);
            x    = x - alpha*mHat./(sqrt(vHat)+obj.epsilon);
        end

    end

    methods (Access = private, Static)

        function value = readOption(s,name,default)
            value = default;
            if isfield(s,name)
                value = s.(name);
            end
        end

    end

end