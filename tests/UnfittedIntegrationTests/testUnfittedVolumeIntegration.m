classdef testUnfittedVolumeIntegration < testUnfittedGeometricalIntegration   
    methods (Access = protected)
        function P = computeGeometricalVariable(obj)
            P = obj.mesh.computeVolume();
        end
    end
end

