classdef testUnfittedGeometricalIntegration < handle
    properties (GetAccess = protected, SetAccess = protected)
        mesh
    end
    
    properties (Access = protected, Abstract)
        topOpt
        meshType
        meshIncludeBoxContour
        levelSet
    end
    
    methods (Access = protected)
        function createMesh(obj)
            mesh_background = obj.topOpt.mesh;
            obj.mesh = Mesh_Unfitted_Factory.create(obj.meshType,mesh_background,Interpolation.create(mesh_background,'LINEAR'),'includeBoxContour',obj.meshIncludeBoxContour);
            obj.mesh.computeMesh(obj.levelSet);
        end
        
        function M = computeGeometricalVariable(obj)
            M = obj.mesh.computeMass();
        end
    end
end

