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

        function C = evaluate(obj,xV)
            lambdaEv = obj.lambda.evaluate(xV);
            muEv     = obj.mu.evaluate(xV);
            nGaus = size(xV,2);
            nElem = length(muEv);
            nStre = 3;
            C = zeros(nStre,nStre,nGaus,nElem);
            C(1,1,:,:)= 2*muEv+lambdaEv;
            C(1,2,:,:)= lambdaEv;
            C(2,1,:,:)= lambdaEv;
            C(2,2,:,:)= 2*muEv+lambdaEv;
            C(3,3,:,:)= muEv;
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.lambda = cParams.lambdaField;
            obj.mu     = cParams.muField;
        end
    end
end