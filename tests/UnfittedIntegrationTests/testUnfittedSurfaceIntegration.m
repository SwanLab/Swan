classdef testUnfittedSurfaceIntegration < testUnfittedGeometricalIntegration
    methods (Access = protected)
        function A = computeGeometricalVariable(obj)
            %             A = obj.mesh.computeSurface();
            A = obj.mesh.computeMass();
        end
    end
end

