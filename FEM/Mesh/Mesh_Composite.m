classdef Mesh_Composite < handle
    
    properties (GetAccess = public, SetAccess = public)
        meshes
        nMeshes
        
        meshInterior
        boxFaceMeshes
        meshBackground
        
        globalConnectivities
        unfittedType
        
        nActiveMeshes
        activeMeshesList
    end
    
    methods (Access = public)
        
        function obj = Mesh_Composite(mesh)
            obj.activeMeshesList = find([false true(1,mesh.nBoxActiveMeshes)]);
            obj.nActiveMeshes = numel(obj.activeMeshesList);
            for iMesh = 1:mesh.nMeshes
                obj.meshes{iMesh} = mesh.meshes{iMesh}.meshBackground;
                obj.meshes{iMesh}.unfittedType = 'INTERIOR';
                obj.meshes{iMesh}.meshBackground = obj.meshes{iMesh};
            end
            obj.meshInterior = mesh.meshInterior;
            obj.globalConnectivities = mesh.globalConnectivities;
            obj.unfittedType = 'COMPOSITE';
            obj.meshBackground = mesh.meshBackground;
            
            obj.createMeshes();
        end
        
        function aMeshes = getActiveMeshes(obj)
            aMeshes = cell(1,obj.nActiveMeshes);
            for iActive = 1:obj.nActiveMeshes
                iMesh = obj.activeMeshesList(iActive);
                aMeshes{iActive} = obj.meshes{iMesh};
            end
        end
    end
    
    methods (Access = private)
        
        function createMeshes(obj)
            
        end
        
    end
    
end