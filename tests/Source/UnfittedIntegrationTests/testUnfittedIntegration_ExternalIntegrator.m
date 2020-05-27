classdef testUnfittedIntegration_ExternalIntegrator < testUnfittedIntegration
    
    methods (Access = protected)
        
        function totalIntegral = computeGeometricalVariable(obj)
            totalIntegral = obj.mesh.computeMass2();
        end
        
    end
    
end

