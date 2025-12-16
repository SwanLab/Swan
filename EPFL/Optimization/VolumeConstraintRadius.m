classdef VolumeConstraintRadius < handle

    properties (Access = private)
        mesh
        volumeTarget
        volume
    end
    
    methods (Access = public)
        function obj = VolumeConstraintRadius(cParams)
            obj.init(cParams);
        end
        
        function [J,dJ] = computeFunctionAndGradient(obj,x)
            [V,dV] = obj.volume.computeFunctionAndGradient(x);
            J      = obj.computeFunction(V);
            dJ{1}  = obj.computeGradient(dV{1});
        end  
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh         = cParams.mesh;
            obj.volumeTarget = cParams.volumeTarget;
            obj.volume       = VolumeFunctionalRadius(cParams);
        end

        function J = computeFunction(obj,V)
            vTar = obj.volumeTarget;
            J    = V/vTar-1;
        end

        function dJ = computeGradient(obj,dV)
            vTar = obj.volumeTarget;
            dJ   = dV;
            dJ.setFValues(dV.fValues/vTar);
        end
    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Volume constraint';
        end
    end
end

