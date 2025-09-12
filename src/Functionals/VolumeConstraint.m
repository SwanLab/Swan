classdef VolumeConstraint < handle

    properties (Access = private)
        mesh
        volumeTarget
        volume
        value0
    end
    
    methods (Access = public)
        function obj = VolumeConstraint(cParams)
            obj.init(cParams);
        end
        
        function [J,dJ] = computeFunctionAndGradient(obj,x)
            x = x{1};
            [V,dV] = obj.volume.computeFunctionAndGradient(x);
            J      = obj.computeFunction(V);
            dJ{1}  = obj.computeGradient(dV);
        end  

        function vol = getDesignVariable(obj)
            vol = obj.volume.getDesignVariable();
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
            if isempty(obj.value0)
              obj.value0 = 1;
            end
            J = J/obj.value0;
        end

        function dJ = computeGradient(obj,dV)
            vTar    = obj.volumeTarget;
            fValues = dV.fValues/vTar;
            fValues = fValues/obj.value0;
            dJ      = FeFunction.create(dV.order,fValues,obj.mesh);
        end
    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Volume constraint';
        end
    end
end

