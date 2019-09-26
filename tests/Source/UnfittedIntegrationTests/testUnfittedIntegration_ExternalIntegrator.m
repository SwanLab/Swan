classdef testUnfittedIntegration_ExternalIntegrator < testUnfittedIntegration
    
    properties (Access = protected)
        integrator
    end
    
    properties (Access = private)
        nodalVector
        integratedVector
    end
    
    methods (Access = protected)
        
        function totalIntegral = computeGeometricalVariable(obj)
            obj.createNodalVector();
            obj.integrateNodalVector();
            totalIntegral = sum(obj.integratedVector);
        end
        
    end
    
    methods (Access = private)
        
        function createNodalVector(obj)
            npnod = obj.mesh.meshBackground.npnod;
            f = ones(npnod,1);
            obj.nodalVector = f;
        end
        
        function int = integrateNodalVector(obj)
            cParams.mesh = obj.mesh;
            cParams.type = 'COMPOSITE';
            switch obj.meshType
                case 'INTERIOR'
                    cParams = obj.createInteriorParams(cParams,obj.mesh);
                case 'BOUNDARY'
                    cParams.compositeParams = cell(0);
                    if contains(class(obj),'Rectangle') || contains(class(obj),'Cylinder')
                        cParams.boxFaceToGlobal = obj.mesh.nodesInBoxFaces;
                        for iMesh = 1:obj.mesh.nActiveBoxFaces
                            iActive = obj.mesh.activeBoxFaceMeshesList(iMesh);
                            boxFaceMesh = obj.mesh.boxFaceMeshes{iActive};
                            params = obj.createCompositeParams(boxFaceMesh);
                            params.boxFaceToGlobal = obj.mesh.nodesInBoxFaces{iActive};
                            cParams.compositeParams{iMesh} = params;
                        end
                    end
                    cParamsInnerCut = obj.createInnerCutParams(obj.mesh);
                    cParams.compositeParams{end+1} = cParamsInnerCut;
            end
            integratorC = Integrator.create(cParams);
            f = obj.nodalVector;
            int = integratorC.integrateAndSum(f);
            obj.integratedVector = int;
        end
        
        function cParams = createInnerCutParams(obj,mesh)
            cParams.mesh = mesh.innerCutMesh;
            cParams.type = 'CutMesh';
        end
        
        function cParams = createInnerParams(obj,mesh)
            cParams.mesh = mesh.innerMesh;
            cParams.type = 'SIMPLE';
            cParams.globalConnec = mesh.globalConnec;
            cParams.npnod = mesh.innerMesh.npnod;
            cParams.backgroundMesh = mesh.meshBackground;
            cParams.innerToBackground = mesh.backgroundFullCells;
        end
        
        function cParams = createCompositeParams(obj,mesh)
            cParams.mesh = mesh;
            cParams.type = 'COMPOSITE';
            cParams.backgroundMesh = obj.mesh.meshBackground;
            cParams.globalConnec = mesh.globalConnec;
            cParams.innerToBackground = [];
            cParams.npnod = obj.mesh.meshBackground.npnod;
            cParams = obj.createInteriorParams(cParams,mesh);
        end
        
        function cParams = createInteriorParams(obj,cParams,mesh)
            cParamsInnerCut = obj.createInnerCutParams(mesh);
            cParams.compositeParams{1} = cParamsInnerCut;
            cParamsInner = obj.createInnerParams(mesh);
            cParams.compositeParams{2} = cParamsInner;
        end
    end
    
end

