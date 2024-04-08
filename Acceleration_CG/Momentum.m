classdef Momentum < handle

    properties (Access = public)
        beta
    end

    properties (Access = private)
        xOld
        momentumCase
        computeBeta
        nIter
    end

    methods (Access = public)

        function obj = Momentum(cParams)
            obj.init(cParams);
        end

        function [y,xF] = apply(obj,x)
            b         = obj.computeBeta();
            y         = x + b*(x - obj.xOld);
            xF        = obj.setEvaluationPoint(x,y);
            obj.xOld  = x;
            obj.nIter = obj.nIter + 1;
            obj.beta  = b;
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.xOld         = cParams.x0;
            obj.momentumCase = cParams.momentumCase;
            obj.nIter        = 1;
            switch cParams.betaStrategy
                case 'Constant'
                    obj.computeBeta = cParams.momentumVal;
                case 'Adaptative'
                    obj.computeBeta = @obj.computeAdaptativeBeta;
                otherwise
                    error('Beta update strategy not implemented')
            end
            obj.beta = obj.computeBeta();
        end
        
        function b = computeAdaptativeBeta(obj)
            k = obj.nIter;
            b = k/(k+2);
        end

        function xF = setEvaluationPoint(obj,x,y)
            switch obj.momentumCase
                case 'Polyak'
                    xF = x;
                case 'Nesterov'
                    xF = y;
                otherwise
                    error('Momentum case not implemented.')
            end
        end

    end

end