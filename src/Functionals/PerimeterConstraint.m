classdef PerimeterConstraint < handle

    properties (Access = private)
        epsilon
        minEpsilon
        target
        target0
        perimeter
        value0
        valueOld
        tarVolume
        totalVolume
    end
    
    methods (Access = public)
        function obj = PerimeterConstraint(cParams)
            obj.init(cParams);
            obj.computeTotalVolume(cParams);
        end
        
        function [J,dJ] = computeFunctionAndGradient(obj,x)
            xD = x.obtainDomainFunction();
            V  = Integrator.compute(xD{1},x.fun.mesh,2);
            V  = V/(obj.tarVolume*obj.totalVolume) - 1;
            obj.updateEpsilonForNextIteration(x.fun,V);
            [P,dP] = obj.perimeter.computeFunctionAndGradient(x);
            J      = obj.computeFunction(P);
            dJ     = obj.computeGradient(dP);
        end  
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.epsilon    = cParams.epsilon;
            obj.minEpsilon = cParams.minEpsilon;
            obj.target     = cParams.target;
            obj.target0    = cParams.target0;
            obj.perimeter  = PerimeterFunctional(cParams);
            obj.value0     = cParams.value0;
            obj.valueOld   = -inf;
            obj.tarVolume  = cParams.tarVolume;
        end

        function computeTotalVolume(obj,s)
            u = ConstantFunction.create(1,s.mesh);
            obj.totalVolume = Integrator.compute(u,s.mesh,2);
        end

        function J = computeFunction(obj,P)
            pTar = obj.target0;
            J    = P/(pTar/obj.value0) - 1; % P-pTar/obj.value0 if pTar is close to zero!!
        end

        function dJ = computeGradient(obj,dP)
            pTar = obj.target0;
            dJ   = dP;
            dJ{1}.setFValues(dP{1}.fValues/(pTar/obj.value0));
        end

        function updateEpsilonForNextIteration(obj,x,V) % Cuando la suma de grays empieza a decaer puede provocar tmb la decay de epsilon
            %if abs(J)<=1e-2
            if norm(x.fValues-obj.valueOld)==0 && V<0.1
                obj.epsilon = obj.epsilon/1.01;
                obj.epsilon = max(obj.epsilon,obj.minEpsilon);
                obj.perimeter.updateEpsilon(obj.epsilon);
                obj.target0 = max(obj.target0*0.8,obj.target);
            end
            obj.valueOld = x.fValues;
            %end % Será preferible tener una decay constante al inicio y luego más notoria hacia el final (cuando el volumen esta por cumplirse y tenemos muchos grises)
        end
    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Perimeter constraint';
        end
    end
end