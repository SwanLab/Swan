classdef PerimeterConstraint < handle

    properties (Access = private)
        epsilon
        minEpsilon
        target
        perimeter
        value0
        valueOld
    end
    
    methods (Access = public)
        function obj = PerimeterConstraint(cParams)
            obj.init(cParams);
        end
        
        function [J,dJ] = computeFunctionAndGradient(obj,x)
            [P,dP] = obj.perimeter.computeFunctionAndGradient(x);
            J      = obj.computeFunction(P);
            dJ     = obj.computeGradient(dP);
            obj.updateEpsilonForNextIteration(J);
        end  
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.epsilon    = cParams.epsilon;
            obj.minEpsilon = cParams.minEpsilon;
            obj.target     = cParams.target;
            obj.perimeter  = PerimeterFunctional(cParams);
            obj.value0     = cParams.value0;
            obj.valueOld   = -inf;
        end

        function J = computeFunction(obj,P)
            pTar = obj.target;
            J    = P/(pTar/obj.value0) - 1; % P-pTar/obj.value0 if pTar is close to zero!!
        end

        function dJ = computeGradient(obj,dP)
            pTar = obj.target;
            dJ   = dP;
            dJ{1}.setFValues(dP{1}.fValues/(pTar/obj.value0));
        end

        function updateEpsilonForNextIteration(obj,J) % Cuando la suma de grays empieza a decaer puede provocar tmb la decay de epsilon
            %if abs(J)<=1e-2
            if J-obj.valueOld<-1e-2 || abs(J)<=1e-2
                obj.epsilon = obj.epsilon/1.01;
                obj.epsilon = max(obj.epsilon,obj.minEpsilon);
                obj.perimeter.updateEpsilon(obj.epsilon);
            end
            obj.valueOld = J;
            %end % Será preferible tener una decay constante al inicio y luego más notoria hacia el final (cuando el volumen esta por cumplirse y tenemos muchos grises)
        end
    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Perimeter constraint';
        end
    end
end