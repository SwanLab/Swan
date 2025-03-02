classdef HeavisideProjector < handle

    properties (Access = private)
        beta
        eta
    end

    methods (Access = public)

        function obj = HeavisideProjector(cParams)
            obj.init(cParams);
        end

        function projF = project(obj,xF)
            a     = tanh(obj.beta*obj.eta) + obj.computeExpressionInNum(xF);
            b     = obj.computeExpressionInDen();
            projF = a/b;
        end

        function derProjF = derive(obj,xF)
            a        = 1 - obj.computeExpressionInNum(xF).^2;
            b        = obj.computeExpressionInDen();
            derProjF = obj.beta*a/b;
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.beta = cParams.beta;
            obj.eta  = cParams.eta;
        end

        function num = computeExpressionInNum(obj,xF)
            num = tanh(obj.beta*(xF.fValues-obj.eta));
        end

        function den = computeExpressionInDen(obj)
            den = tanh(obj.beta*obj.eta) + tanh(obj.beta*(1-obj.eta));
        end
        
    end    
end