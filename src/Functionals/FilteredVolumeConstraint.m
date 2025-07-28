classdef FilteredVolumeConstraint < handle

    properties (Access = private)
        mesh
        alpha
        volume
    end
    
    methods (Access = public)
        function obj = FilteredVolumeConstraint(cParams)
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
            obj.mesh   = cParams.mesh;
            obj.alpha  = cParams.alpha;
            obj.volume = FilteredVolumeFunctional(cParams);
        end

        function J = computeFunction(obj,V)
            volFrac = obj.alpha;
            J       = V/volFrac-1;
        end

        function dJ = computeGradient(obj,dV)
            volFrac = obj.alpha;
            dJ      = dV./volFrac;
        end
    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Filtered volume constraint';
        end
    end
end

