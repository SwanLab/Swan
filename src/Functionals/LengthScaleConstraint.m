classdef LengthScaleConstraint < handle

    properties (Access = private)
        volume
        perimeter
        target
        target0
        epsilon
        minEpsilon
        delta
        valueOld
    end

    methods (Access = public)
        function obj = LengthScaleConstraint(cParams)
            obj.init(cParams);
        end

        function [J,dJ] = computeFunctionAndGradient(obj,x)
            [V,dV] = obj.volume.computeFunctionAndGradient(x);
            V      = V*obj.volume.totalVolume;
            [P,dP] = obj.perimeter.computeFunctionAndGradient(x);
            h0     = obj.target0;
            h      = obj.target;
            t      = obj.delta;
            J      = -(V+2*t*h)/((P+t)*h0) + 1;
            dJ{1}  = copy(dV{1});
            dV     = dV{1}.fValues*obj.volume.totalVolume;
            dP     = dP{1}.fValues;
            dJ{1}.setFValues(-dV./((P+t)*h0) + dP.*((V+2*t*h)/((P+t)^2*h0)));
            obj.updateForNextIteration(J);
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.volume     = VolumeFunctional(cParams);
            obj.perimeter  = PerimeterFunctional(cParams);
            obj.target     = cParams.target;
            obj.target0    = obj.target*100;
            obj.epsilon    = cParams.epsilon;
            obj.minEpsilon = cParams.minEpsilon;
            obj.delta      = cParams.delta;
            obj.valueOld   = -inf;
        end

        function updateForNextIteration(obj,J)
            if J-obj.valueOld<0 || abs(J)<=1e-2
                obj.epsilon = obj.epsilon/1.01;
                obj.epsilon = max(obj.epsilon,obj.minEpsilon);
                obj.perimeter.updateEpsilon(obj.epsilon);
                obj.target0 = max(obj.target0/1.008,obj.target);
            end
            obj.valueOld = J;
        end
    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'LS Constr';
        end
    end
end