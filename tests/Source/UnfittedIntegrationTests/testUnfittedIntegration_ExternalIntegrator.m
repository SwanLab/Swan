classdef testUnfittedIntegration_ExternalIntegrator < testUnfittedIntegration
    
    methods (Access = protected)
        
        function totalIntegral = computeGeometricalVariable(obj)
            switch obj.meshType  
                case 'INTERIOR'
                    totalIntegral = obj.unfittedMesh.computeMass();
                case 'BOUNDARY'
                    totalIntegral = obj.unfittedMesh.computePerimeter();
            end

        end
        
    end
    
end

