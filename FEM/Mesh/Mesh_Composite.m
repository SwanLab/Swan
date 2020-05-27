classdef Mesh_Composite < AbstractMesh
    
    properties (GetAccess = public, SetAccess = public)
        meshes
        nMeshes
        nActiveMeshes
        activeMeshesList
        
        globalConnectivities
    end
    
    methods (Access = public)
        
        function obj = Mesh_Composite()
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
    
    methods (Access = protected)
        
        function append(obj,mesh)
           obj.meshes{end+1} = mesh;
        end
        
    end
    
    methods (Access = private)
        
        function createMeshes(obj)
            
        end
        
    end
    
end