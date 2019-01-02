classdef testUnfittedGeometricalIntegration < handle
    properties (GetAccess = protected, SetAccess = protected)
        mesh
    end
    
    properties (Access = protected, Abstract)
        topOpt
        meshType
        levelSet
    end
    
    methods (Access = protected)
        function createMesh(obj)
            mesh_background = obj.topOpt.mesh;
            obj.mesh = Mesh_Unfitted(obj.meshType,mesh_background,Interpolation.create(mesh_background,'LINEAR'));
            obj.mesh.computeMesh(obj.levelSet);
        end
    end
end

