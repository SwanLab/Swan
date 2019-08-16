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
        
        function int = integrateMesh(obj)
            %             int = obj.integrateMeshOld();
            int = obj.integrateMeshNew();
        end
        
        function M2 = integrateMeshOld(obj)
            cParams.mesh = obj.mesh;
            cParams.type = cParams.mesh.unfittedType;
            obj.integrator = Integrator.create(cParams);
            M2 = obj.integrator.integrate(ones(size(obj.mesh.levelSet_background)));
        end
        
        function int = integrateMeshNew(obj)
            switch obj.meshType
                case 'INTERIOR'
                    cParams.mesh = obj.mesh;
                    cParams.type = 'COMPOSITE';
                    cParamsInnerCut = obj.createInnerCutParams(obj.mesh);
                    cParams.compositeParams{1} = cParamsInnerCut;
                    cParamsInner = obj.createInnerParams(obj.mesh);
                    cParams.compositeParams{2} = cParamsInner;
                case 'BOUNDARY'
                    if contains(class(obj),'Rectangle') || contains(class(obj),'Cylinder')
                        cParams.type = 'COMPOSITE';
                        cParams.mesh = obj.mesh;
                        cParams.boxFaceToGlobal = obj.mesh.nodesInBoxFaces;
                        for iMesh = 1:obj.mesh.nActiveBoxFaces
                            iActive = obj.mesh.activeBoxFaceMeshesList(iMesh);
                            boxFaceMesh = obj.mesh.boxFaceMeshes{iActive};
                            params = obj.createCompositeParams(boxFaceMesh);
                            params.boxFaceToGlobal =  obj.mesh.nodesInBoxFaces{iActive};
                            cParams.compositeParams{iMesh} = params;
                        end
                        cParamsInnerCut = obj.createInnerCutParams(obj.mesh);
                        cParams.compositeParams{end+1} = cParamsInnerCut;
                    else
                        cParams.mesh = obj.mesh;
                        cParams.type = 'COMPOSITE';
                        cParamsInnerCut = obj.createInnerCutParams(obj.mesh);
                        cParams.compositeParams{1} = cParamsInnerCut;
                    end
            end
            integratorC = Integrator.create(cParams);
            f = ones(size(obj.mesh.levelSet_background));
%             int = integratorC.integrate(f);
            int = integratorC.integrateAndSum(f);
        end
        
        function params = createInnerCutParams(obj,mesh)
            params.mesh = mesh.innerCutMesh;
            params.type = 'CutMesh';
        end
        
        function params = createInnerParams(obj,mesh)
            params.mesh = mesh.innerMesh;
            params.type = 'SIMPLE';
            params.globalConnec = mesh.globalConnec;
            params.npnod = mesh.innerMesh.npnod;
            params.backgroundMesh = mesh.meshBackground;
            params.innerToBackground = mesh.backgroundFullCells;
        end
        
        function params = createCompositeParams(obj,mesh)
            params.mesh = mesh;
            params.type = 'COMPOSITE';
            params.backgroundMesh = obj.mesh.meshBackground;
            params.globalConnec = mesh.globalConnec;
            params.innerToBackground = [];
            params.npnod = obj.mesh.meshBackground.npnod;
            cParamsInner = obj.createInnerParams(mesh);
            params.compositeParams{1} = cParamsInner;
            cParamsInnerCut = obj.createInnerCutParams(mesh);
            params.compositeParams{2} = cParamsInnerCut;
        end
    end
    
end

