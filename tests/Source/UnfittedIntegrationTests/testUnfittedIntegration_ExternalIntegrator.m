classdef testUnfittedIntegration_ExternalIntegrator < testUnfittedIntegration
    
    properties (Access = protected)
        integrator
    end
    
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
            cParams.mesh = obj.mesh;
            obj.integrator = Integrator.create(cParams);
            %M2 = obj.integrator.integrateUnfittedMesh(ones(size(obj.mesh.levelSet_background)),obj.mesh);
            M2 = obj.integrator.computeIntegral(ones(size(obj.mesh.levelSet_background)));
        end
        
    end
    
end

