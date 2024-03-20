classdef VolumeConstraint < handle

    properties (Access = private)
        mesh
        volumeTarget
        volume
    end
    
    methods (Access = public)
        function obj = VolumeConstraint(cParams)
            obj.init(cParams);
        end
        
        function [J,dJ] = computeFunctionAndGradient(obj,x)
            [V,dV] = obj.volume.computeFunctionAndGradient(x);
            J      = obj.computeFunction(V);
            dJ     = obj.computeGradient(dV);
        end  
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh         = cParams.mesh;
            obj.volumeTarget = cParams.volumeTarget;
            obj.volume       = VolumeFunctional(cParams);
        end

        function J = computeFunction(obj,V)
            vTar = obj.volumeTarget;
            J    = V/vTar-1;
        end

        function dJ = computeGradient(obj,dV)
            vTar       = obj.volumeTarget;
            for iF = 1:numel(dV)
                dJ{iF}         = dV{iF}.copy();
                dJ{iF}.fValues = dV{iF}.fValues/vTar;
            end
        end
    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Volume constraint';
        end
    end
end

