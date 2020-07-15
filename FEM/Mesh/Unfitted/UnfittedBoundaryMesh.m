classdef UnfittedBoundaryMesh < handle
    
    properties (Access = public)
        meshes
        globalConnec
    end
    
    properties (Access = private)
        nBoundaries
        levelSet
        activeMeshes
    end
    
    properties (Access = private)
        boundaryMesh
    end
    
    methods (Access = public)
        
        function obj = UnfittedBoundaryMesh(cParams)
            obj.init(cParams)
            obj.createUnfittedMeshes();
            obj.obtainGlobalConnec()
        end
        
        function compute(obj,ls)
            obj.levelSet = ls;
            obj.computeActiveMesh();
            obj.computeUnfittedMeshes();
        end
        
        function m = getActiveMesh(obj)
            m = obj.getActiveField(obj.meshes);                         
        end
        
        function g = getGlobalConnec(obj)
            g = obj.getActiveField(obj.globalConnec);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.boundaryMesh = cParams.boundaryMesh;
            obj.nBoundaries  = numel(obj.boundaryMesh);
            obj.activeMeshes = false(obj.nBoundaries,1);            
        end
        
        function createUnfittedMeshes(obj)
            for iBoundary = 1:obj.nBoundaries
                s.backgroundMesh = obj.boundaryMesh{iBoundary}.mesh;
                s.isInBoundary   = true;
                s = SettingsMeshUnfitted(s);
                obj.meshes{iBoundary} = UnfittedMesh(s);                
            end
        end
        
        function obtainGlobalConnec(obj)
            for iBoundary = 1:obj.nBoundaries
                m = obj.boundaryMesh{iBoundary};
                obj.globalConnec{iBoundary} = m.nodesInBoxFaces;
            end            
        end
        
        function computeActiveMesh(obj)
            for iBoundary = 1:obj.nBoundaries                
                isActive = obj.isUnfittedMeshActive(iBoundary);
                obj.activeMeshes(iBoundary) = isActive;
            end            
        end
        
        function computeUnfittedMeshes(obj)
            for iBoundary = 1:obj.nBoundaries       
                isMeshActive = obj.activeMeshes(iBoundary);
                if isMeshActive
                   obj.computeUnfittedMesh(iBoundary);
                end
            end            
        end     
        
        function computeUnfittedMesh(obj,iBoundary)
            nodes = obj.globalConnec{iBoundary};
            ls    = obj.levelSet(nodes);            
            obj.meshes{iBoundary}.compute(ls);
        end
        
        function itIs = isUnfittedMeshActive(obj,iBoundary)
            nodes = obj.globalConnec{iBoundary};
            ls    = obj.levelSet(nodes);          
            itIs = any(sign(ls)<0);
        end
        
        function activeFields = getActiveField(obj,fields)
            activeFields = cell(0);
            iMesh = 1;
            for iBoundary = 1:obj.nBoundaries       
                isMeshActive = obj.activeMeshes(iBoundary);
                if isMeshActive
                   activeFields{iMesh} = fields{iBoundary};
                   iMesh = iMesh +1;
                end
            end              
        end        
        
    end
    
end