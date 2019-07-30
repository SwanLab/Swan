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
                        cParams.type = cParams.mesh.unfittedType;
            obj.integrator = Integrator.create(cParams);
            %M2 = obj.integrator.integrateUnfittedMesh(ones(size(obj.mesh.levelSet_background)),obj.mesh);
            M2 = obj.integrator.integrate(ones(size(obj.mesh.levelSet_background)));
             
%             cParams.mesh = obj.mesh;
%             cParams.type = 'COMPOSITE';
%             interiorIntegrator = Integrator.create(cParams);
%             M2_2 = interiorIntegrator.integrate(ones(size(obj.mesh.levelSet_background)));
%             sum(abs(M2-M2_2{1}))
%             
        end
        
    end
    
end

