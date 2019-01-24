classdef testUnfittedIntegration_ExternalIntegrator < testUnfittedIntegration
    
    methods (Access = protected)
        
        function M = computeGeometricalVariable(obj)
            M2 = obj.integrateMesh();
            M = obj.sumResults(M2);
        end
        
        function B = sumResults(~,A)
            B = sum(A);
        end
        
    end
    
    methods (Access = private)
        
        function M2 = integrateMesh(obj)
            integrator = Integrator.create(obj.mesh);
            M2 = integrator.integrateUnfittedMesh(ones(size(obj.mesh.x_background)),obj.mesh);
        end
        
    end
    
end

