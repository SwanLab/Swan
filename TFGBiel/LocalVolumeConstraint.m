classdef LocalVolumeConstraint < handle

    properties (Access = private)
        epsilon
        minEpsilon
        target
        volume
        value0
    end
    
    methods (Access = public)
        function obj = LocalVolumeConstraint(cParams)
            obj.init(cParams);
        end
        
        function [J,dJ] = computeFunctionAndGradient(obj,x)
            [P,dP]  = obj.volume.computeFunctionAndGradient(x);
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
            obj.volume     = LocalVolumeFunctional(cParams);
            obj.value0     = cParams.value0;
        end

        function J = computeFunction(obj,P)
            pTar = obj.target;
            J    = P-pTar/obj.value0;
        end

        function dJ = computeGradient(obj,dj)
%             pTar = obj.target;
            dj.setFValues(dj.fValues);
            dJ{1} = dj;
        end

        function updateEpsilonForNextIteration(obj,J)
            if abs(J)<=1e-2
                obj.epsilon = obj.epsilon/1.001;
                obj.epsilon = max(obj.epsilon,obj.minEpsilon);
                obj.volume.updateEpsilon(obj.epsilon);
            end
        end
    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Volume constraint';
        end
    end
end