classdef testUnfittedVolumeIntegration < testUnfittedGeometricalIntegration
    methods (Access = protected)
        function V = computeGeometricalVariable(obj)
            %             V = obj.mesh.computeVolume();
            V = obj.mesh.computeMass();
        end
    end
end

