classdef InnerMesh < handle
    
    properties (GetAccess = public, SetAccess = private)
       fullCells
       backgroundMesh
       mesh
       globalConnec
    end
    
    properties (Access = private)
        all2unique
        unique2all
        uniqueNodes
        coord
        connec
    end
        
    methods (Access = public)
        
        function obj = InnerMesh(cParams)
            obj.init(cParams);
            obj.computeUniqueNodes();
            obj.computeCoords();
            obj.computeConnec();
            obj.createMesh();
            obj.computeGlobalConnec();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.backgroundMesh = cParams.backgroundMesh;
            obj.fullCells      = cParams.fullCells;
        end
        
        function computeUniqueNodes(obj)
            connecB = obj.backgroundMesh.connec;
            connecF = connecB(obj.fullCells,:);
            allNodes = connecF(:);
            [uNodes,ind,ind2] = unique(allNodes,'rows','stable');
            obj.all2unique  = ind;
            obj.unique2all  = ind2;
            obj.uniqueNodes = uNodes;
        end
        
        function computeCoords(obj)
            uNodes       = obj.uniqueNodes;
            allCoords    = obj.backgroundMesh.coord;
            uniqueCoords = allCoords(uNodes,:);
            obj.coord    = uniqueCoords;
        end
        
        function computeConnec(obj)
            nnode = size(obj.backgroundMesh.connec,2);
            nCell = size(obj.fullCells,1);
            obj.connec = reshape(obj.unique2all,nCell,nnode);
        end
        
        function createMesh(obj)
            s.coord  = obj.coord;
            s.connec = obj.connec;
            s.kFace  = obj.backgroundMesh.kFace;
            obj.mesh = Mesh.create(s);
        end
        
        function computeGlobalConnec(obj)
            con  = obj.backgroundMesh.connec;
            obj.globalConnec = con(obj.fullCells,:);
        end
        
    end
    
end