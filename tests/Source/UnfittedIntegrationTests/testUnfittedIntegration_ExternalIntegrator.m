classdef testUnfittedIntegration_ExternalIntegrator < testUnfittedIntegration
    
    methods (Access = protected)
        
        function totalIntegral = computeGeometricalVariable(obj)
            switch obj.meshType  
                case 'INTERIOR'
                    totalIntegral = obj.mesh.computeMass();
                case 'BOUNDARY'
                    totalIntegral = obj.mesh.computePerimeter();
            end

        end
        
    end
    
end

