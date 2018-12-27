classdef testUnfittedPerimeterIntegration < testUnfittedGeometricalIntegration   
    methods (Access = protected)
        function P = computeGeometricalVariable(obj)
            P = obj.mesh.computePerimeter();
        end
    end
end

