classdef testUnfittedPerimeterIntegration < testUnfittedGeometricalIntegration
    methods (Access = protected)
        function P = computeGeometricalVariable(obj)
            %             P = obj.mesh.computePerimeter();
            P = obj.mesh.computeMass();
        end
    end
end

