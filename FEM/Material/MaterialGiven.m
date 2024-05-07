classdef MaterialGiven < handle

    properties (Access = private)
        lambda
        mu
        matValues
    end

    methods (Access = public)
        function obj = MaterialGiven(cParams)
            obj.init(cParams);
        end

        function matEv = evaluate(obj,xV)
            lambdaEv = obj.lambda.evaluate(xV);
            muEv     = 
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.lambda = cParams.lambdaField;
            obj.mu     = cParams.muField;
        end
    end
end