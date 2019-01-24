classdef testUnfittedIntegration_InternalIntegrator < testUnfittedIntegration
    
    methods (Access = protected)
        
        function M = computeGeometricalVariable(obj)
            M = obj.mesh.computeMass();
        end
        
    end
    
end

