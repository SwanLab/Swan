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
%             cParamsInnerCut = obj.createInnerCutParams();
%             cParams.compositeParams{1} = cParamsInnerCut;
%             cParamsInner = obj.createInnerParams();
%             cParams.compositeParams{2} = cParamsInner;
%             integratorC = Integrator.create(cParams);
%             M2_2 = integratorC.integrate(ones(size(obj.mesh.levelSet_background)));
%             ref = sum(M2);
%             cut = sum(M2_2{1});
%             inner = sum(M2_2{2});
%             total = M2_2{1} + M2_2{2};
%             sum(abs(M2-total))
%             
        end
        
        function params = createInnerCutParams(obj)
            params.mesh = obj.mesh.innerCutMesh;
            params.type = 'CutMesh';
        end
        
        function params = createInnerParams(obj)
            params.mesh = obj.mesh.innerMesh;
            params.type = 'SIMPLE';
            params.globalConnec = obj.mesh.globalConnec;
            params.npnod = obj.mesh.innerMesh.npnod;
            params.backgroundMesh = obj.mesh.meshBackground;
            params.innerToBackground = obj.mesh.backgroundFullCells;
        end
    end
    
end

