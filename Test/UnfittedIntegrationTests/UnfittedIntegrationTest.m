classdef UnfittedIntegrationTest < handle
    
    properties (Access = private)
        mesh
        meshType
        levelSet
        analyticalValue
    end

    properties (Access = private)
        unfittedMesh
        varAdim
    end
    
    methods (Access = public)
        
        function obj = UnfittedIntegrationTest(cParams)
            obj.init(cParams);
            obj.integrateSurface();
        end
        
        function error = computeError(obj)
            error = abs(obj.varAdim - 1);
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.mesh            = cParams.mesh;
            obj.meshType        = cParams.meshType;
            obj.levelSet        = cParams.levelSet;
            obj.analyticalValue = cParams.analyticalValue;
        end
        
        function integrateSurface(obj)
            obj.createMesh();
            geomVar     = obj.computeGeometricalVariable();
            obj.varAdim = geomVar/obj.analyticalValue;
        end

        function createMesh(obj)
            s.backgroundMesh = obj.mesh;
            s.boundaryMesh   = obj.mesh.createBoundaryMesh();
            obj.unfittedMesh = UnfittedMesh(s);
            obj.unfittedMesh.compute(obj.levelSet);
        end

        function totalIntegral = computeGeometricalVariable(obj)
            switch obj.meshType
                case 'INTERIOR'
                    totalIntegral = obj.unfittedMesh.computeMass();
                case 'BOUNDARY'
                    totalIntegral = obj.unfittedMesh.computePerimeter();
                case 'SUBBOUNDARY' % Provisional
                    totalIntegral = 0;
                    meshes        = obj.unfittedMesh.unfittedBoundaryMesh.getActiveMesh();
                    nBoundary     = size(meshes,2);
                    uniqueCoord   = cell(0);
                    for i = 1:nBoundary
                        subUnfMesh    = meshes{i};
                        subMeshes     = subUnfMesh.unfittedBoundaryMesh.getActiveMesh();
                        nSubMeshes    = size(subMeshes,2);
                        for j = 1:nSubMeshes
                            coordj        = subMeshes{j}.backgroundMesh.coord;
                            prevCoord     = cell2mat(uniqueCoord);
                            if i>1
                                isRepeated    = all(ismember(coordj,prevCoord,'rows'));
                                if(isRepeated)
                                    subMeshes{j}.compute(ones(1000,1));
                                else
                                    uniqueCoord   = [uniqueCoord;coordj];
                                end
                            else
                                uniqueCoord   = [uniqueCoord;coordj];
                            end
                        end
                        totalIntegral = totalIntegral + subUnfMesh.computePerimeter();
                    end
            end
        end

    end

end