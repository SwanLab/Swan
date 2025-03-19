classdef PerimeterConstraint < handle

    properties (Access = private)
        epsilon
        minEpsilon
        target
        perimeter
    end
    
    methods (Access = public)
        function obj = PerimeterConstraint(cParams)
            obj.init(cParams);
        end
        
        function [J,dJ] = computeFunctionAndGradient(obj,x)
            [P,dP]  = obj.perimeter.computeFunctionAndGradient(x);
            J       = obj.computeFunction(P);
            dJ      = obj.computeGradient(dP{1});
            obj.updateEpsilonForNextIteration(J);
        end  
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.epsilon    = cParams.epsilon;
            obj.minEpsilon = cParams.minEpsilon;
            obj.target     = cParams.target;
            obj.perimeter  = PerimeterFunctional(cParams);
        end

        function J = computeFunction(obj,P)
            pTar = obj.target;
            J    = P/pTar-1;
        end

        function dJ = computeGradient(obj,dJ)
            pTar = obj.target;
            dJ.setFValues(dJ.fValues./pTar);
        end

        function updateEpsilonForNextIteration(obj,J)
            if abs(J)<=1e-2
                obj.epsilon = obj.epsilon/1.001;
                obj.epsilon = max(obj.epsilon,obj.minEpsilon);
                obj.perimeter.updateEpsilon(obj.epsilon);
            end
        end
    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Perimeter constraint';
        end
    end
end