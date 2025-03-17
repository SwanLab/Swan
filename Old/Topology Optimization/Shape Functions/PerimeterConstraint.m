classdef PerimeterConstraint < handle

    properties (Access = private)
        mesh
        epsilon
        minEpsilon
        perimeterTargetAbs
        perimeter
    end
    
    methods (Access = public)
        function obj = PerimeterConstraint(cParams)
            obj.init(cParams);
        end
        
        function [J,dJ] = computeFunctionAndGradient(obj,x)
            [P,dP]  = obj.perimeter.computeFunctionAndGradient(x);
            J       = obj.computeFunction(P);
            dJ      = obj.computeGradient(dP);
            obj.updateEpsilonForNextIteration(J);
        end  
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh               = cParams.mesh;
            obj.epsilon            = cParams.epsilon;
            obj.minEpsilon         = cParams.minEpsilon;
            obj.perimeterTargetAbs = cParams.perimeterTargetAbs;
            cParams.value0         = 1;
            obj.perimeter          = PerimeterFunctional(cParams);
        end

        function J = computeFunction(obj,P)
            pTar = obj.perimeterTargetAbs;
            J    = P/pTar-1;
        end

        function dJ = computeGradient(obj,dP)
            pTar    = obj.perimeterTargetAbs;
            fValues = dP.fValues/pTar;
            dJ      = FeFunction.create(dP.order,fValues,obj.mesh);
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