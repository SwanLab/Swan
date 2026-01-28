classdef VolumeConstraintParameter < handle

    properties (Access = private)
        mesh
        volumeTarget
        volume
    end
    
    methods (Access = public)
        function obj = VolumeConstraintParameter(cParams)
            obj.init(cParams);
        end
        
        function [J,dJ] = computeFunctionAndGradient(obj,x)
            [V,dV] = obj.volume.computeFunctionAndGradient(x);
            J      = obj.computeFunction(V);
            dJ  = obj.computeGradient(dV);
        end  
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh         = cParams.mesh;
            obj.volumeTarget = cParams.volumeTarget;
            obj.volume       = VolumeFunctionalParameter(cParams);
        end

        function J = computeFunction(obj,V)
            vTar = obj.volumeTarget;
            J    = V/vTar-1;
        end

        function dJ = computeGradient(obj,dV)
            nf = length(dV);
            vTar = obj.volumeTarget;
            for i = 1:nf
                dJ{i}   = dV{i};
                dJ{i}.setFValues(dJ{i}.fValues/vTar);
            end
        end
    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Volume constraint';
        end
    end
end

