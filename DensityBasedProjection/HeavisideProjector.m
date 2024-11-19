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
            if isempty(obj.eta)
                projF = 1 - exp(-obj.beta*xF.fValues) + xF.fValues*exp(-obj.beta);
            else
                a     = tanh(obj.beta*obj.eta) + obj.computeExpressionInNum(xF);
                b     = obj.computeExpressionInDen();
                projF = a/b;
            end
        end

        function derProjF = derive(obj,xF)
            if isempty(obj.eta)
                derProjF = obj.beta*exp(-obj.beta*xF.fValues) + exp(-obj.beta);
            else
                a        = 1 - obj.computeExpressionInNum(xF).^2;
                b        = obj.computeExpressionInDen();
                derProjF = obj.beta*a/b;
            end
        end


        function updateBeta(obj, beta)
            obj.beta = beta;
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.beta = cParams.beta;
            if isfield(cParams, 'eta')
                obj.eta  = cParams.eta;
            end
        end

        function num = computeExpressionInNum(obj,xF)
            num = tanh(obj.beta*(xF.fValues-obj.eta));
        end

        function den = computeExpressionInDen(obj)
            den = tanh(obj.beta*obj.eta) + tanh(obj.beta*(1-obj.eta));
        end

    end    
end